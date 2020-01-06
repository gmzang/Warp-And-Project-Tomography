#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include<iostream>
#include<algorithm>

#include<scope.h>

template< typename real >
real ls_abs( const real &in ){
    return (in < real(0.0) ? -in : in);
}

template< typename real, typename real3, typename distance_func >
real3 ls_signed_distance_gradient( distance_func &f, const real3 &query, const real eps=1e-3 ){
    const real f0 = f( query );
    const real f1 = f( query+real3(eps,0,0) );
    const real f2 = f( query+real3(0,eps,0) );
    const real f3 = f( query+real3(0,0,eps) );
    return real3( f1-f0, f2-f0, f3-f0 )/eps;
}

template< typename real, typename real3, typename F1, typename F2 >
inline real ls_signed_distance_union( const F1& f1, const F2& f2, const real3& query ){
	real v1=f1(query);
	real v2=f2(query);
	return std::min( v1, v2 );
}

template< typename real, typename real3, typename F1, typename F2 >
inline real ls_signed_distance_intersection( const F1& f1, const F2& f2, const real3& query ){
	return std::max( f1(query), f2(query) );
}

template< typename real, typename real3, typename F1, typename F2 >
inline real ls_signed_distance_difference( const F1& f1, const F2& f2, const real3& query ){
	return std::max( f1(query), -f2(query) );
}

template< typename real, typename real3, typename F >
inline real ls_signed_distance_complement( const F& f, const real3& query ){
	return -f(query);
}

template< typename real, typename real3 >
real ls_signed_distance_point( const real &point, const real3 &query ){
    return query-point.length();
}

// negative inside the sphere
template< typename real, typename real3 >
real ls_signed_distance_sphere( const real3 &center, const real r, const real3 &query ){
    return (query-center).length()-r;
}

// negative inside the box
template< typename real, typename real3 >
real ls_signed_distance_box( const real3 &minim, const real3 &maxim, const real3 &query ){
    real3 clamp = query;
    if( query[0] >= minim[0] && query[0] <= maxim[0] && query[1] >= minim[1] && query[1] <= maxim[1] && query[2] >= minim[2] && query[2] <= maxim[2] ){
        // inside the box
        clamp[0] = std::min( ls_abs(query[0]-minim[0]), ls_abs(query[0]-maxim[0]) );
        clamp[1] = std::min( ls_abs(query[1]-minim[1]), ls_abs(query[1]-maxim[1]) );
        clamp[2] = std::min( ls_abs(query[2]-minim[2]), ls_abs(query[2]-maxim[2]) );
        return - std::min( clamp[0], std::min( clamp[1], clamp[2] ) );
    }
    // outside the box, can clamp the query point to the
    // bounding box to get the closest point (since the
    // box is orthogonal).
    clamp[0] = std::max( minim[0], std::min( maxim[0], query[0] ) );
    clamp[1] = std::max( minim[1], std::min( maxim[1], query[1] ) );
    clamp[2] = std::max( minim[2], std::min( maxim[2], query[2] ) );
    return (query-clamp).length();
}

/**
 @brief expression-templated level set base class, all this class serve to do is
 act as a base class for defining level set operator in terms of the algebraic
 operators. without it the compiler can't tell that it's supposed to be operating
 on level sets.  it provides casting operators to the underlying level-set type
 which in turn inherit from this class via the curiously-recursive template 
 design pattern.
 
 in order to determine types, all level set classes must inherit from this class,
 define the real and real3 data types as well as providing a call operator taking
 a real3 point in world space and returning the signed distance from the level
 set as a result
*/
template< typename LS >
class ls_expr {
public:
	inline operator const LS&() const {
		return static_cast<const LS&>(*this);
	}
	inline operator LS&(){
		return static_cast<const LS&>(*this);
	}
};

/**
 @brief defines a sphere as a subclass of the level-set expression class
 ls_expr above.
*/
template< typename real_type, typename real3_type>
class ls_sphere : public ls_expr< ls_sphere<real_type,real3_type> > {
public:
	typedef real_type  real;
	typedef real3_type real3;
private:
	real3  m_center;
	real   m_radius;
public:
	ls_sphere( const real3 center, const real radius ) : m_center(center), m_radius(radius) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_sphere<real,real3>( m_center, m_radius, p );
	}
};

/**
 @brief defines a box as a subclass of the level-set expression class
 ls_expr above.
*/
template< typename real_type, typename real3_type >
class ls_box : public ls_expr< ls_box<real_type,real3_type> > {
public:
	typedef real_type  real;
	typedef real3_type real3;
private:
	real3 m_minim;
	real3 m_maxim;
public:
	ls_box( const real3 minim, const real3 maxim ) : m_minim(minim), m_maxim(maxim) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_box<real,real3>( m_minim, m_maxim, p );
	}
};

/**
 @brief defines an image as a subclass of the level-set expression class
 ls_expr above
*/
template< typename image_type, typename real_type, typename real3_type >
class ls_image : public ls_expr< ls_image<image_type,real_type,real3_type> > {
public:
	typedef real_type real;
	typedef real3_type real3;
	typedef image_type image;
private:
	int		m_dim[3];
	real	m_aabb[6];
	image	m_img;
public:
	ls_image( const int *dim, const real *aabb, const image& img ){
		for( int i=0; i<3; i++ ) m_dim[i] = dim[i];
		for( int i=0; i<6; i++ ) m_aabb[i] = aabb[i];
		m_img = img;
	}
	inline real operator()( const real3 &p ) const {
		real3 g( real(m_dim[0])*(p[0]-m_aabb[0])/(m_aabb[1]-m_aabb[0]),
				 real(m_dim[1])*(p[1]-m_aabb[2])/(m_aabb[3]-m_aabb[2]),
				 real(m_dim[2])*(p[2]-m_aabb[4])/(m_aabb[5]-m_aabb[4]) );
		return m_img.linear_atXYZ( g[0], g[1], g[2] );
	}
};

/**
 @brief defines the level-set union, this simply calls the level-set API function
 ls_signed_distance_union using the two templated types, note that the constructor
 to this function takes ls_expr references, which are then down-cast to the 
 actual type via the cast operator in ls_expr
*/
template< typename LS1, typename LS2 >
class ls_union : public ls_expr< ls_union< LS1, LS2 > > {
public:
	typedef typename LS1::real  real;
	typedef typename LS2::real3 real3;
private:
	const LS1& m_ls1;
	const LS2& m_ls2;
public:
	ls_union( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ) : m_ls1(ls1), m_ls2(ls2) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_union<real,real3,LS1,LS2>( m_ls1, m_ls2, p );
	}
};

template< typename LS1, typename LS2 >
ls_union< LS1, LS2 > operator+( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ){
	return ls_union<LS1,LS2>( ls1, ls2 );
}

/**
 @brief defines the level-set difference, this simply calls the level-set API function
 ls_signed_distance_difference using the two templated types, note that the constructor
 to this function takes ls_expr references, which are then down-cast to the
 actual type via the cast operator in ls_expr
 */
template< typename LS1, typename LS2 >
class ls_difference : public ls_expr< ls_difference< LS1, LS2 > > {
public:
	typedef typename LS1::real  real;
	typedef typename LS2::real3 real3;
private:
	const LS1& m_ls1;
	const LS2& m_ls2;
public:
	ls_difference( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ) : m_ls1(ls1), m_ls2(ls2) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_difference<real,real3,LS1,LS2>( m_ls1, m_ls2, p );
	}
};
template< typename LS1, typename LS2 >
ls_difference< LS1, LS2 > operator-( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ){
	return ls_difference<LS1,LS2>( ls1, ls2 );
}

/**
 @brief defines the level-set intersection, this simply calls the level-set API function
 ls_signed_distance_intersection using the two templated types, note that the constructor
 to this function takes ls_expr references, which are then down-cast to the
 actual type via the cast operator in ls_expr
 */
template< typename LS1, typename LS2 >
class ls_intersection : public ls_expr< ls_difference< LS1, LS2 > > {
public:
	typedef typename LS1::real  real;
	typedef typename LS2::real3 real3;
private:
	const LS1& m_ls1;
	const LS2& m_ls2;
public:
	ls_intersection( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ) : m_ls1(ls1), m_ls2(ls2) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_intersection<real,real3,LS1,LS2>( m_ls1, m_ls2, p );
	}
};
template< typename LS1, typename LS2 >
ls_difference< LS1, LS2 > operator*( const ls_expr<LS1>& ls1, const ls_expr<LS2>& ls2 ){
	return ls_intersection<LS1,LS2>( ls1, ls2 );
}

/**
 @brief defines the level-set complement, this simply calls the level-set API function
 ls_signed_distance_complement using the two templated types, note that the constructor
 to this function takes ls_expr references, which are then down-cast to the
 actual type via the cast operator in ls_expr
 */
template< typename LS >
class ls_complement : public ls_expr< ls_complement< LS > > {
public:
	typedef typename LS::real  real;
	typedef typename LS::real3 real3;
private:
	const LS& m_ls;
public:
	ls_complement( const ls_expr<LS>& ls ) : m_ls(ls) {
	}
	inline real operator()( const real3& p ) const {
		return ls_signed_distance_complement<real,real3,LS>( m_ls, p );
	}
};

template< typename LS >
ls_complement<LS> operator-( const ls_expr<LS>& ls ){
	return ls_complement<LS>( ls );
}

#endif
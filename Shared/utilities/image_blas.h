#ifndef IMAGE_BLAS_H
#define IMAGE_BLAS_H

#include<cmath>
#include<cassert>

namespace image_blas {

	template< typename real, typename image >
	class scaled_image {
	private:
		const image &m_image;
		const real   m_scale;
	public:
		scaled_image( const image &img, const real scale ) : m_image(img), m_scale(scale) { }
		inline int width() const { return m_image.width(); }
		inline int height() const { return m_image.height(); }
		inline int depth() const { return m_image.depth(); }
		inline int spectrum() const { return m_image.spectrum(); }
		inline real operator()( const int i, const int j, const int k, const int c ) const {
			return m_image(i,j,k,c)*m_scale;
		}
	};
	
	template< typename real,typename image >
	scaled_image<real,image> scaled( const image &img, const real scale ){
		return scaled_image<real,image>(img,scale);
	}
	
	template< typename I >
	void copy( const I &a, I &b ){
		b=a;
	}
	
	template< typename real, typename I1, typename I2 >
	real dot( const I1 &a, const I2 &b ){
		const int dim[] = { a.width(), a.height(), a.depth(), a.spectrum() };
		assert( b.width()==dim[0] && b.height()==dim[1] && b.depth()==dim[2] && b.spectrum()==dim[3] );
		real res = real(0.0);
		for( int i=0; i<dim[0]; i++ ){
			for( int j=0; j<dim[1]; j++ ){
				for( int k=0; k<dim[2]; k++ ){
					for( int c=0; c<dim[3]; c++ ){
						res += a(i,j,k,c)*b(i,j,k,c);
					}
				}
			}
		}
		return res;
	}
	
	template< typename real, typename I >
	real norm( const I &a ){
		return sqrt(dot( a, a ));
	}
	
	template< typename real, typename I >
	real norm_squared( const I &a ){
		return dot(a,a);
	}
	
	template< typename I1, typename I2 >
	void add( const I1 &a, I2 &b ){
		const int dim[] = { a.width(), a.height(), a.depth(), a.spectrum() };
		assert( b.width()==dim[0] && b.height()==dim[1] && b.depth()==dim[2] && b.spectrum()==dim[3] );
		for( int i=0; i<dim[0]; i++ ){
			for( int j=0; j<dim[1]; j++ ){
				for( int k=0; k<dim[2]; k++ ){
					for( int c=0; c<dim[3]; c++ ){
						b(i,j,k,c) += a(i,j,k,c);
					}
				}
			}
		}
	}
	
	template< typename I1, typename I2, typename I3 >
	void add( const I1 &a, const I2 &b, I3 &c ){
		const int dim[] = { a.width(), a.height(), a.depth(), a.spectrum() };
		assert( b.width()==dim[0] && b.height()==dim[1] && b.depth()==dim[2] && b.spectrum()==dim[3] );
		c = I3(dim[0],dim[1],dim[2],dim[3]);
		for( int i=0; i<dim[0]; i++ ){
			for( int j=0; j<dim[1]; j++ ){
				for( int k=0; k<dim[2]; k++ ){
					for( int m=0; m<dim[3]; m++ ){
						c(i,j,k,m) = b(i,j,k,m) + a(i,j,k,m);
					}
				}
			}
		}
	}
	
};

#endif
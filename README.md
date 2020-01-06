# Warp-And-Project-Tomography
Source code for SIGGRAPH 2019 Paper: Warp-and-Project Tomography for Rapidly Deforming Objects


[**Warp-And-Project-Tomography**](https://vccimaging.org/Publications/Zang2019WarpAndProject/)

<img src=./teaser.png width=95%/>
<img src=./bubble.png width=95%/>
<img src=./smoke.png width=95%/>
We introduce a CT reconstruction method for objects that undergo rapid deformation during the scan. Shown here is a copper foam crumpling under a compressive force during the scan. The whole complex animation is reconstructed using only 192 projection images that all correspond to different deformation states of the foam.

## Usage
The code is tested in Visual Studio 2015 and 2018 on 64 bits Windows 7 and Windows 10.

To run the Warp-and-Project tomography framwork (WaPTomo), first compile `WaP19.sln` with Visual Studio. Notice that libraries like [argtable], [openmp],  [eigen], and [cimg] are required, which are all included in `./shared` folder, put them in the right path for the project.

More details for running parameters, just type:
`WaP19  --help`

## Software to release

We are developing a open source software for 3D and 4D tomograhic reconstruction, in which the super resolution tomography, space-time tomography, and warp-and-project will be intergrated with a decent GUI (Qt-based front-end and C++ based back-end), more memory/time friendly code. The version 1.0 is currently under finalization and will be released very soon, stay tuned!

<img src=./software.png width=95%/>





## Contact: 
If you find any bug or if you have any suggestion or comment, please contact: 

Guangming : guangming.zang@kaust.edu.sa

Copyright: [Visual computing center], CEMSE, KAUST


[openmp]: <http://openmp.org/wp/>
[eigen]: <http://eigen.tuxfamily.org/index.php?title=Main_Page>
[argtable]: <http://argtable.sourceforge.net/>
[cimg]: <http://cimg.eu/>
[Visual computing center]: <https://vcc.kaust.edu.sa/Pages/Home.aspx>



## License and citation
This research is released under the [CC BY-NC 3.0 US license](https://creativecommons.org/licenses/by-nc/3.0/us/). We encourage an attitude of reproducible work for academic-only purpose. Please kindly cite our work in your publications if it helps your research:

```

@article{zang2019warp,
  title={Warp-and-Project Tomography for Rapidly Deforming Objects},
  author={Zang, Guangming and Idoughi, Ramzi and Tao, Ran and Lubineau, Gilles and Wonka, Peter and Heidrich, Wolfgang},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={4},
  year={2019},
  publisher={ACM} 
  } 

```

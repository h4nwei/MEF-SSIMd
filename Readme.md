
# Perceptual Evaluation for Multi-Exposure Image Fusion of Dynamic Scenes
This is the implementation for [Perceptual Evaluation for Multi-Exposure Image Fusion of Dynamic Scenes](sim.jxufe.cn/JDMKL/pdf/19_TIP_MEF-SSIMd.pdf), Yuming Fang, Hanwei Zhu, Kede Ma, Zhou Wang, Shutao Li, IEEE Transactions on Image Processing (TIP), to appear, 2019.

## Abstract
A common approach to high dynamic range (HDR) imaging is to capture multiple images of different exposures followed by multi-exposure image fusion (MEF) in either radiance or intensity domain. A predominant problem of this approach is
the introduction of the ghosting artifacts in dynamic scenes with camera and object motion. While many MEF methods (often referred to as deghosting algorithms) have been proposed for reduced ghosting artifacts and improved visual quality, little work has been dedicated to perceptual evaluation of their deghosting
results. Here we first construct a database that contains 20 multi-exposure sequences of dynamic scenes and their corresponding fused images by nine MEF algorithms. We then carry out a subjective experiment to evaluate fused image quality, and find that none of existing objective quality models for MEF provides
accurate quality predictions. Motivated by this, we develop an objective quality model for MEF of dynamic scenes. Specifically, we divide the test image into static and dynamic regions, measure structural similarity between the image and the corresponding sequence in the two regions separately, and combine quality
measurements of the two regions into an overall quality score. Experimental results show that the proposed method significantly outperforms the state-of-the-art. In addition, we demonstrate the promise of the proposed model in parameter tuning of MEF methods.

## Prerequisites

- MATLAB

## Dataset

The MEF image quality database can be obtained at: http://sim.jxufe.cn/JDMKL/code/DeghostingIQADatabase.rar. The dataset contains:

- 20 multi-exposure image sequences of dynamic natural scenes
- 9 MEF methods  generate 180 fused images
- Mean opinion scores of the fused images


## Test

Here, we test the source image sequnce 'horse' with two fused images, and show the quality maps.
 ```
run demo.m
 ```

 ## Reference
 - Y. Fang, H. Zhu, K. Ma, and Z. Wang, “Perceptual quality assessment of HDR deghosting algorithms,” in *IEEE International Conference on Image Processing*, 2017, pp. 3165–3169.
  - K. Ma, K. Zeng, and Z. Wang, “Perceptual quality assessment for multi-exposure image fusion,” *IEEE Transactions on Image Processing*, vol. 24, no. 11, pp. 3345–3356, Nov. 2015.

 ## Citation
 ```
 @article{Fang2019,
  title={Perceptual Evaluation for Multi-Exposure Image Fusion of Dynamic Scenes},
  author={Yuming Fang, Hanwei Zhu, Kede Ma, Zhou Wang, Shutao Li},
  journal={IEEE Transactions on Image Processing},
  publisher={to appear}
}

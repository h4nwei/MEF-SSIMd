# Perceptual Evaluation of Multi-Exposure Image Fusion for Dynamic Scene
This is the implementation for [Perceptual Evaluation of Multi-Exposure Image Fusion for Dynamic Scene](sim.jxufe.cn/JDMKL/code/DeghostingIQADatabase.rar), Yuming Fang, Hanwei Zhu, Kede Ma, Zhou Wang, Shutao Li, IEEE Transactions on Image Processing (TIP), to appear, 2019.

## Abstract
A common approach to high dynamic range (HDR) imaging is to capture multiple images of different exposures followed by multi-exposure image fusion (MEF) in either radiance or intensity domain. A predominant problem of this approach is
the introduction of the ghosting artifacts in dynamic scenes with camera and object motion. While many MEF methods (often referred to as deghosting algorithms) have been proposed for reduced ghosting artifacts and improved visual quality, little work has been dedicated to perceptual evaluation of their deghosting
results. Here we first construct a database that contains 20 multi-exposure sequences of dynamic scenes and their corresponding fused images by nine MEF algorithms. We then carry out a subjective experiment to evaluate fused image quality, and find that none of existing objective quality models for MEF provides
accurate quality predictions. Motivated by this, we develop an objective quality model for MEF of dynamic scenes. Specifically, we divide the test image into static and dynamic regions, measure structural similarity between the image and the corresponding sequence in the two regions separately, and combine quality
measurements of the two regions into an overall quality score. Experimental results show that the proposed method significantly outperforms the state-of-the-art. In addition, we demonstrate the promise of the proposed model in parameter tuning of MEF methods.

## Prerequisites
----
- MATLAB

## Dataset
----
The contruction MEF image quality database can be obtained at: sim.jxufe.cn/JDMKL/code/DeghositngIQADatabase.rar. The dataset contains
```

```



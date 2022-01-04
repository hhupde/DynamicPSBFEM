Free and Forced Vibration Analysis in Abaqus Based on the Polygonal Scaled Boundary Finite Element Method


The polygonal scaled boundary finite element method (PSBFEM) is a novel approach integrating the standard scaled boundary finite element method and the polygonal mesh technique. In this work, a user-defined element (UEL) for dynamic analysis based on the PSBFEM is developed in the general finite element software ABAQUS. We present the main procedures of interacting with Abaqus, updating AMATRX and RHS, defining the UEL element, and solving the stiffness and mass matrices through eigenvalue decomposition. Several benchmark problems of free and forced vibration are solved to validate the proposed implementation. The results show that the PSBFEM is more accurate than the FEM with mesh refinement. Moreover, the PSBFEM avoids the occurrence of hanging nodes by constructing a polygonal mesh. Thus, the PSBFEM can choose an appropriate mesh resolution for different structures ensuring accuracy and reducing calculation costs.

Citation
If you use our code for academic research, you are encouraged to cite the following paper:

```
@article{ye2021free,
  title={Free and Forced Vibration Analysis in Abaqus Based on the Polygonal Scaled Boundary Finite Element Method},
  author={Ye, Nan and Su, Chao and Yang, Yang},
  journal={Advances in Civil Engineering},
  volume={2021},
  year={2021},
  publisher={Hindawi}
}

```

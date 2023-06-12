
# Fin analysis using finite element method (FEM)
- github repository for computational physics, final project [2023.06.14]

# how to use it
1. [import FEM_2D](#1.-import-FEM-2D)

2. [set mesh's boundary](#2.-set-mesh's-boundary)

3. [build mesh](#3-build-mesh)

4. [set mesh's geometry](#4-set-meshs-geometry)

5. [set boundary condition](#5-set-boundary-condition)

6. [compile the model](#6-compile-the-model)

7. [compute mesh](#7-compute-mesh)

8. [save and import result](#8-save-and-import-result)

## 1. import FEM 2D

FEM_2D folder should be exists in current directory like below.

<p align="center">
    <img src="./images/import%20FEM_2d.png" style="height :150px" alt="Current file is main.ipynb" title="imoprt FEM 2D"/>
</p>

In main.ipynb, you can import FEM_2D like below.

```python
from FEM_2D import *

Test_fin = FEM_2D(*args)
```

## 2. set mesh's boundary

Before you define mesh, you should define the boundary.

```python
Test_fin = FEM_2(X=x_max, Y=y_max, dx=unit_x, dy=unit_y)
```

Here `X` represent the maximum boundary of mesh along X axis, and `Y` represent the maximum of Y axis.

`dx` and `dy` represent the unit length of mesh, so the size of one element of mesh will be `dx * dy`. It would be better to equalize `dx`, `dy` in most case.

So when you set `X=10, Y=15, dx=0.1, dy=0.1`, your mesh will be set like this.

<p align="center">
    <img src="./images/empty_mesh.png" style="height :300px" alt="X=10, Y=15, dx=0.1, dy=0.1" title="X=10, Y=15, dx=0.1, dy=0.1"/>
</p>



## 3. build mesh

## 4. set mesh's geometry

## 5. set boundary condition

## 6. compile the model

## 7. compute mesh

## 8. save and import result



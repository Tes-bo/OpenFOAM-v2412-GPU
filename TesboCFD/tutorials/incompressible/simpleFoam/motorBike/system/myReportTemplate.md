---
marp: true
paginate: true
---

<style>
:root {
    font-size: 20px;
}
td {
    width: 1000px;
}
table {
    width: 100%;
}
img {
    display: block;
    margin-left: auto;
    margin-right: auto;
    width: 60%;
}
</style>

# {{OF_EXECUTABLE}} : {{OF_CASE_NAME}} tutorial

- Case: {{OF_CASE_PATH}}
- Submission: {{OF_CLOCK_START}} on {{OF_DATE_START}}
- Report time: {{OF_CLOCK_NOW}} on {{OF_DATE_NOW}}

---

## Run information

| Property       | Value              |
|----------------|--------------------|
| Host           | {{OF_HOST}}        |
| Processors     | {{OF_NPROCS}}      |
| Time steps     | {{OF_TIME_INDEX}}  |
| Initial deltaT | {{initial_deltaT}} |
| Current deltaT | {{OF_TIME_DELTAT}} |
| Execution time | {{executionTime}}  |

---

## OpenFOAM information

| Property       | Value              |
|----------------|--------------------|
| Version        | {{OF_VERSION}}     |
| API            | {{OF_API}}         |
| Patch          | {{OF_PATCH}}       |
| Build          | {{OF_BUILD}}       |
| Architecture   | {{OF_BUILD_ARCH}}  |

---

## Mesh statistics

| Property          | Value                |
|-------------------|----------------------|
| Bounds            | {{OF_MESH_BOUNDS_MIN}}{{OF_MESH_BOUNDS_MAX}} |
| Number of cells   | {{OF_MESH_NCELLS}}   |
| Number of faces   | {{OF_MESH_NFACES}}   |
| Number of points  | {{OF_MESH_NPOINTS}}  |
| Number of patches | {{OF_MESH_NPATCHES}} |

---

## Linear solvers

| Property | Value          | tolerance(rel)   | Tolerance(abs)      |
|----------|----------------|------------------|---------------------|
| p        | `{{solver_p}}` | {{solver_p_tol}} | {{solver_p_reltol}} |
| U        | `{{solver_U}}` | {{solver_u_tol}} | {{solver_u_reltol}} |

---

## Numerical scehemes

The chosen divergence schemes comprised:

~~~
{{divSchemes}}
~~~

---

## Graphs

Residuals

![]({{OF_CASE_PATH}}/postProcessing/residualGraph1/{{OF_TIME}}/residualGraph1.svg)

---

## Results

Forces

![]({{OF_CASE_PATH}}/postProcessing/forceCoeffsGraph1/{{OF_TIME}}/forceCoeffsGraph1.svg)

---

Made using Open&nabla;FOAM v2412 from https://openfoam.com

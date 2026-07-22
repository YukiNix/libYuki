# libYuki

libYuki is a personal scientific coding library of Yuki Nix / 禰夕・有希 的个人日用科研代码库. 

## Code Layout / 代码结构

```text
julia/
    libYukiBasic.jl
    libYukiConstant.jl
    libYukiMath.jl
        libYukiMathVector.jl
        libYukiMathStatistics.jl
        libYukiMathInterpolation.jl
    libYukiPhysics.jl
        libYukiPhysicsBody.jl
        libYukiPhysicsField.jl
        libYukiPhysicsNBodySimulation.jl
    libYukiAstronomy.jl
        libYukiAstronomyTransit.jl
            libYukiAstronomyTransitLimbDarkening.jl
            libYukiAstronomyTransitTTV.jl
        libYukiAstronomyKeplerMission.jl
```

## Main Capabilities / 主要功能

- Basic helpers and measurement utilities / 基础工具与测量值处理
- Constants (physics, astrophysics, mathematics) / 常数库(物理, 天体物理, 数学)
- Numerical math (vector, interpolation, statistics, differentiation, integration) / 数值数学(向量, 插值, 统计, 微分, 积分)
- Physics and simulation helpers (including N-body workflows) / 物理计算与模拟辅助(含 N 体流程)
- Domain-specific toolchains used in daily work, including transit and Kepler-related utilities / 日常工作中使用的领域化工具链(含凌星与 Kepler 相关工具)

## Dependencies / 依赖

Current Julia code uses / 当前 Julia 代码使用:

- CSV
- DataFrames
- Dates
- DifferentialEquations
- FITSIO
- ForwardDiff
- HTTP
- Interpolations
- JLD2
- LinearAlgebra
- Measurements
- MultivariateStats
- Optim
- Orbits
- ProgressMeter
- PyCall
- QuadGK
- Random
- StableRNGs
- StaticArrays
- Statistics
- Tables
- Transits
- Turing

Install example / 安装示例:

```julia
using Pkg
Pkg.add([
        "CSV", "DataFrames", "DifferentialEquations", "FITSIO", "ForwardDiff",
        "HTTP", "Interpolations", "JLD2", "Measurements", "MultivariateStats",
        "Optim", "Orbits", "ProgressMeter", "PyCall", "QuadGK", "StableRNGs",
        "StaticArrays", "Tables", "Transits", "Turing"
])
```

## Loading / 加载方式

This repository is include-based (not a registered Julia package yet) / 当前以 include 脚本组织(尚未注册为 Julia 包):

```julia
include("julia/libYukiAstronomy.jl")
```

## Notes / 说明

- APIs are still evolving while modules are being consolidated / API 仍在持续整合中. 
- Some files include each other directly; loading through entry files is recommended / 部分文件会互相 include，建议通过入口文件加载. 
- Example notebooks in repository root are no longer actively updated / 根目录示例 notebook 已停止更新. 

## License / 许可证

BSD 3-Clause. See LICENSE / BSD 3-Clause，详见 LICENSE. 


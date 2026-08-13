# nav-adjust-go-lib

`nav-adjust-go-lib` 是一个独立的 Go 站网平差库。它接收站点之间互相观测得到的相对 ENU 向量及可选的控制点坐标先验，通过加权最小二乘计算统一基准下的站点 ENU 坐标、平差后基线、残差与质量指标。

该库不依赖 RTK 库，也不接收 RTK 内部结构。调用方只需按本库定义的 `Station`、`ENUBaseline` 和 `ENUNetworkProblem` 组织数据。

## 核心模型

一条从站点 `i` 指向站点 `j` 的观测定义为：

```text
ObservedENU(i→j) = PositionENU(j) - PositionENU(i) + error
```

每条观测同时包含 `East`、`North`、`Up` 三个分量及其完整 `3×3` 协方差矩阵。库求解：

```text
min vᵀ R⁻¹ v
```

其中 `v = Observed - Adjusted`，`R` 是所有基线和坐标先验协方差组成的随机模型。

必须遵守以下数据约定：

- 所有站点坐标和相对向量都采用同一套 ENU 轴方向、同一原点定义和同一单位，推荐使用米；
- `From=A, To=B` 表示向量 `B - A`，反向使用时三个分量都必须取反；
- `Covariance` 填方差和协方差，不是标准差，也不是权值，矩阵必须有限、对称、正定；
- 默认情况下，每个互不连通的站网分量至少包含一个固定站或控制点坐标先验，否则整体平移无法确定；纯相对站网可以显式启用零质心自由网基准；
- 当前固定站坐标按无误差基准处理。自由站不需要初始坐标。

如果每条基线的 ENU 分量来自各自不同的局部坐标系，必须先在业务层把它们转换到共同坐标系，不能直接送入本库。

## 安装

```bash
go get github.com/wfu-work/nav-adjust-go-lib
```

大型项目推荐显式使用数据模型包和站网求解包：

```go
import (
    "github.com/wfu-work/nav-adjust-go-lib/model"
    "github.com/wfu-work/nav-adjust-go-lib/network"
)
```

简单项目也可以只导入根包 `github.com/wfu-work/nav-adjust-go-lib`。根包是兼容门面，类型均为别名，不产生转换和额外分配；两种方式都不会向调用方暴露 Gonum 类型。

## 分包设计

```text
nav-adjust-go-lib (根包门面)
├── model/       公共输入、输出、矩阵和 JSON 数据契约；不依赖求解器
├── network/     ENU 站网校验、拓扑、方程组装、求解和结果映射
├── core/        求解器无关的线性方程、约束和底层结果模型
├── batch/       批处理加权最小二乘及协方差块求解
├── internal/sparse/ CSR 对称稀疏矩阵、投影、Jacobi/块 Jacobi-PCG 求解
├── quality/     卡方检验和残差质量工具
├── nonlinear/   通用 Gauss-Newton 迭代能力
├── robust/      通用 Huber 迭代重加权能力
├── variance/    通用分组协方差尺度更新能力
└── sequential/  通用线性 Kalman 能力
```

依赖方向保持单向：

```text
调用方 → model
调用方 → network → model
                 → batch → core
                 → quality → core
根包门面 → model + network
```

- 只负责存储、传输或 JSON 解析的代码可以仅依赖 `model`；
- 站间 ENU 向量业务直接依赖 `model + network`；
- 根包保留 `SolveENUNetwork` 等完整名称，适合简单调用和向后兼容；
- `core/`、`batch/` 等目录是可复用数值层，不依赖任何 RTK 或站点业务类型。

## 最小示例

```go
package main

import (
    "fmt"

    "github.com/wfu-work/nav-adjust-go-lib/model"
    "github.com/wfu-work/nav-adjust-go-lib/network"
)

func main() {
    datum := model.ENU{East: 0, North: 0, Up: 0}
    covariance := model.Matrix3FromStdDev(0.01, 0.01, 0.02)

    problem := model.Problem{
        Name: "demo-network",
        Stations: []model.Station{
            {ID: "A", Name: "datum", Fixed: true, KnownENU: &datum},
            {ID: "B"},
            {ID: "C"},
        },
        Baselines: []model.Baseline{
            {
                ID: "AB", From: "A", To: "B",
                Vector: model.ENU{East: 10.002, North: 0.001, Up: 0.003},
                Covariance: covariance,
            },
            {
                ID: "AC", From: "A", To: "C",
                Vector: model.ENU{East: 0.001, North: 9.998, Up: -0.002},
                Covariance: covariance,
            },
            {
                ID: "BC", From: "B", To: "C",
                Vector: model.ENU{East: -10.000, North: 10.001, Up: -0.004},
                Covariance: covariance,
            },
        },
    }

    if err := network.Validate(problem, nil); err != nil {
        panic(err)
    }
    result, err := network.Solve(problem, nil)
    if err != nil {
        panic(err)
    }

    for _, station := range result.Stations {
        fmt.Printf("%s: E=%.4f N=%.4f U=%.4f\n",
            station.ID,
            station.Position.East,
            station.Position.North,
            station.Position.Up,
        )
    }
    fmt.Printf("sigma0=%.4f, dof=%d\n",
        result.Diagnostics.Sigma0,
        result.Diagnostics.DegreesOfFreedom,
    )
}
```

如果控制点坐标本身有不确定度，不要把它设为 `Fixed`，应使用带完整 `3×3` 协方差的坐标先验：

```go
problem.Priors = []model.Prior{
    {
        ID:        "control-a",
        StationID: "A",
        Position:  model.ENU{East: 100, North: 200, Up: 3},
        Covariance: model.Matrix3FromStdDev(0.01, 0.01, 0.02),
    },
}
```

`Fixed` 表示完全无误差的精确坐标；`Prior` 表示参与平差、允许产生残差的随机控制坐标。同一个站不能同时为 `Fixed` 并配置 `Prior`。

## 自由网内部基准

只有站间相对 ENU 向量、没有固定站和控制点坐标先验时，站网存在 E/N/U 三个整体平移自由度。默认的 `DatumExternal` 会拒绝这种输入；如果业务只需要内部相对坐标，可以显式采用零质心基准：

```go
result, err := network.Solve(problem, &model.Options{
    Datum: model.DatumFreeCentroid,
})
```

库会对每个没有外部基准的连通分量分别增加三条精确内部约束：

```text
sum(E_i) = 0
sum(N_i) = 0
sum(U_i) = 0
```

已有固定站或坐标先验的分量保持原有外部基准，不会被重新平移。输出的基线、闭合差和残差不受内部原点选择影响，但自由分量的站点坐标及协方差只在该零质心内部基准下有意义，不能直接解释成外部绝对 ENU 坐标。结果会通过 `Diagnostics.DatumMode`、`FreeDatumComponentCount`、`InternalDatumConstraintCount` 和 `Warnings` 明确报告这一点。

零质心基准使用精确约束。请求 `SolverPCG` 时，库通过约束零空间投影直接求解稀疏法方程，`Diagnostics.Solver` 为 `sparse-projected-pcg`，不会为自由网构造稠密 KKT。完全没有基线或坐标先验的孤立站不会因自由网模式而被伪造为已知坐标，仍返回 `ErrDisconnectedNetwork`。

仓库中有可直接运行的版本：

```bash
go run ./examples/enu-network
```

## 输入 API

### Station

```go
type Station struct {
    ID       string
    Name     string
    Fixed    bool
    KnownENU *ENU
    Metadata map[string]string
}
```

- `ID` 在站网内必须唯一且非空；
- 固定站必须提供有限的 `KnownENU`；
- 自由站只需提供 `ID`，其坐标由平差求出；
- `Name` 和 `Metadata` 不参与计算，会原样带入输出。

### ENUBaseline

```go
type ENUBaseline struct {
    ID         string
    From       string
    To         string
    Vector     ENU
    Covariance Matrix3
    Group      string
    Metadata   map[string]string
}
```

对角协方差可使用辅助方法构造：

```go
// 参数是 E/N/U 的标准差，内部自动平方成方差。
covariance := model.Matrix3FromStdDev(0.01, 0.01, 0.02)

// 参数直接是 E/N/U 方差。
covariance = model.DiagonalMatrix3(0.0001, 0.0001, 0.0004)
```

存在分量相关性时使用行主序完整矩阵：

```go
covariance := model.Matrix3{Data: [9]float64{
    varianceE, covarianceEN, covarianceEU,
    covarianceEN, varianceN, covarianceNU,
    covarianceEU, covarianceNU, varianceU,
}}
```

不同基线之间目前按相互独立处理；一条基线内部的 E/N/U 相关性会完整参与计算。

### PositionPrior

```go
type PositionPrior struct {
    ID         string
    StationID  string
    Position   ENU
    Covariance Matrix3
    Metadata   map[string]string
}
```

- `StationID` 必须指向一个自由站；
- 一个自由站可以配置多个相互独立的坐标先验；
- `Covariance` 的规则与基线相同，支持 E/N/U 完整相关性；
- 坐标先验可单独或与精确固定站混合建立站网基准；
- 鲁棒基线降权不会改变控制点先验的权。

## 对外方法

```go
func network.Validate(
    problem model.Problem,
    options *model.Options,
) error

func network.Solve(
    problem model.Problem,
    options *model.Options,
) (*model.Result, error)

func network.SolveContext(
    ctx context.Context,
    problem model.Problem,
    options *model.Options,
) (*model.Result, error)
```

根包兼容入口为 `ValidateENUNetwork`、`SolveENUNetwork` 和 `SolveENUNetworkContext`。`nil` options 表示普通加权最小二乘和默认质量参数。求解器没有全局可变状态；不同调用之间彼此独立。

带 `Context` 的方法会在稀疏迭代、鲁棒重加权、方差分量迭代和按需协方差查询期间协作式响应取消。稠密 Gonum 分解及调用方自己的 `nonlinear.Model.Linearize` 单次调用不能从外部强制中断，但库会在这些调用前后检查取消状态：

```go
ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
defer cancel()

result, err := network.SolveContext(ctx, problem, options)
if errors.Is(err, context.DeadlineExceeded) {
    // 本次平差已超时停止
}
```

## 大型站网与稀疏求解

默认配置仍使用稠密 Cholesky 和完整协方差，保持已有调用行为不变。参数较多时，可显式启用 CSR 稀疏法方程和 PCG。ENU 网络默认使用连续站点 E/N/U 三参数组成的 `3×3` 块 Jacobi 预条件器；也可以显式选择标量 Jacobi：

```go
result, err := network.Solve(problem, &model.Options{
    Solver: model.SolverOptions{
        Method:            model.SolverPCG,
        Preconditioner:    model.PreconditionerBlockJacobi, // ENU 网络默认值
        MaxIterations:     5000,  // 0 表示按参数规模自动选择
        RelativeTolerance: 1e-10, // 0 表示使用默认值
        AbsoluteTolerance: 1e-12, // 0 表示使用默认值
    },
    Covariance: model.CovarianceStationBlocks,
})
```

协方差模式决定内存与诊断完整度：

| 模式 | 顶层 N×N 协方差 | 站点 3×3 协方差 | 基线/先验残差诊断 | 适用场景 |
| --- | --- | --- | --- | --- |
| `CovarianceFull` | 有 | 有 | 有 | 默认兼容模式，中小站网 |
| `CovarianceStationBlocks` | 无 | 有 | 有 | 需要站点精度与观测诊断的大型站网 |
| `CovarianceNone` | 无 | 无 | 无 | 只需要坐标、原始残差和整体统计的超大站网或鲁棒中间迭代 |

`station-blocks` 不构造完整逆矩阵。它通过重复 PCG 求解，仅提取每个站点以及实际基线/先验涉及的逆法方程元素；内存随稀疏法方程和站网边数增长，但计算时间会高于 `none`。`none` 仍返回坐标、原始残差、`Objective`、`Sigma0`、自由度和全局卡方检验；协方差派生的标准差、标准化残差与冗余度保持零值，并由 `Diagnostics.*Available` 明确标识为不可用。

PCG 的 `Diagnostics.ConditionNumberAvailable` 为 `false`，因为迭代法不虚构精确条件数；`SolverPreconditioner` 报告实际使用的预条件器，应结合 `SolverIterations` 和 `SolverRelativeResidual` 判断收敛质量。未达到容差会返回可通过 `errors.Is(err, network.ErrNotConverged)` 判断的错误。存在精确约束时，库只对规模通常很小的约束 Gram 矩阵 `CCᵀ` 做稠密 Cholesky，并在约束零空间内执行稀疏 PCG；坐标求解和按需受约束协方差都不构造参数规模的稠密 KKT。约束线性相关时返回 `ErrRankDeficient`。

底层 `batch.Solve` 同样支持 `batch.Options{Solver: batch.SolverPCG, Covariance: batch.CovarianceNone}`，通用批处理默认使用标量 Jacobi；若参数天然连续成组，可设置 `PreconditionerBlockJacobi` 和正的 `PreconditionerBlockSize`。需要少量局部逆矩阵元素时，可使用 `batch.SolveDetailed` 的 `FormalCovarianceBlock` 或 `FormalCovarianceValues`，而不必生成完整逆矩阵；这些 API 也分别提供 `SolveDetailedContext`、`FormalCovarianceBlockContext` 和 `FormalCovarianceValuesContext`。

## 方差分量估计

当不同设备、解算策略或观测时段给出的协方差尺度不一致时，可以通过 `Baseline.Group` 将基线分组，并让库迭代估计每组协方差的倍率：

```go
problem.Baselines[0].Group = "receiver-a"
problem.Baselines[1].Group = "receiver-b"

result, err := network.Solve(problem, &model.Options{
    VarianceComponents: &model.VarianceComponentOptions{
        MaxIterations: 20,   // 默认 10
        Tolerance:     1e-3, // 相邻两次尺度的最大相对变化
        MinScale:      1e-4,
        MaxScale:      1e4,
    },
    Solver: model.SolverOptions{Method: model.SolverPCG},
    Covariance: model.CovarianceStationBlocks,
})
```

每个不同的 `Group` 是一个方差分量；空字符串也是合法的默认组。迭代使用组内加权残差平方和与有效冗余度更新倍率：

```text
newScale = oldScale × groupObjective / groupRedundancy
```

最终实际使用的基线随机模型为：

```text
EffectiveCovariance = InputCovariance × VarianceScale / RobustWeight
```

结果中的 `VarianceComponents` 按组名排序，给出 `Scale`、`StdDevScale=sqrt(Scale)`、基线数、观测分量数、目标函数贡献和有效冗余度；每条 `BaselineResult.VarianceScale` 也会记录其组倍率。倍率大于 1 表示输入协方差偏乐观，需要放大；小于 1 表示输入协方差偏保守。

每组必须有足够的独立冗余度。树状站网、仅有一条无多余观测的组等情况无法独立估计，会返回 `ErrInsufficientRedundancy`。控制点坐标先验不参与分组尺度估计，始终保留调用方提供的协方差。启用方差分量估计后，组尺度本身就是由当前数据估计的，因此不再返回普通固定随机模型的全局卡方检验；只启用鲁棒估计时，则仅在基线实际发生降权后省略该检验。

估计中间轮次不会构造完整 N×N 逆矩阵；只查询站点和实际观测涉及的协方差元素，所以可与稀疏 PCG 组合。到达迭代上限仍会返回最终结果，并通过 `Diagnostics.VarianceComponentsConverged=false` 和 `Warnings` 明确提示。

## 鲁棒平差

开启 Huber 迭代重加权：

```go
result, err := network.Solve(problem, &model.Options{
    Robust: &model.RobustOptions{
        Method:        model.RobustHuber,
        Threshold:     2.5,
        MaxIterations: 20,
        Tolerance:     1e-3,
        MinWeight:     0.05,
    },
})
```

鲁棒权以整条 ENU 基线为一个块计算，而不是分别任意修改 E/N/U 三个权。这样可以保留该基线协方差矩阵内部的相关结构。

若鲁棒迭代实际对基线降权，权值已由观测残差自适应确定，结果不再满足普通固定权卡方检验的前提；此时 `Diagnostics.GlobalTest` 为 `nil`，并在 `Warnings` 中说明。若所有基线权值仍为 1，则仍返回普通全局卡方检验。`Sigma0` 可作为最终模型的描述性指标。

输出中的 `BaselineResult.Weight < 1` 和 `Downweighted=true` 表示该基线被降权。鲁棒估计用于降低粗差影响，不等于自动证明某条观测错误；是否剔除仍应由业务规则决定。

## 输出结果

`ENUNetworkResult` 主要包含：

| 字段 | 含义 |
| --- | --- |
| `Stations` | 固定站和自由站的最终 ENU 坐标、标准差及 `3×3` 协方差 |
| `Baselines` | 原始向量、平差后向量、残差、标准化残差、冗余度和鲁棒权 |
| `Priors` | 控制点先验坐标、平差坐标、残差、标准化残差和冗余度 |
| `ParameterKeys` | 完整参数协方差矩阵的行列顺序 |
| `FormalCovariance` | 先验单位权方差为 1 时的自由站参数协方差 |
| `Covariance` | `Sigma0² × FormalCovariance`；无多余观测时 `Sigma0=1` |
| `Diagnostics.Objective` | 加权残差平方和 `vᵀR⁻¹v` |
| `Diagnostics.DegreesOfFreedom` | 自由度 |
| `Diagnostics.Sigma0` | 后验单位权中误差 |
| `Diagnostics.SolverIterations` | PCG 参数求解迭代次数；直接法为 0 |
| `Diagnostics.SolverRelativeResidual` | PCG 参数求解最终相对残差 |
| `Diagnostics.SolverPreconditioner` | PCG 实际使用的 `jacobi` 或 `block-jacobi`；直接法为空 |
| `Diagnostics.DatumMode` | 外部基准或零质心自由网基准 |
| `Diagnostics.FreeDatumComponentCount` | 使用零质心内部基准的连通分量数 |
| `Diagnostics.InternalDatumConstraintCount` | 为自由分量加入的内部基准约束数 |
| `VarianceComponents` | 各基线组估计出的协方差倍率、标准差倍率、目标函数与冗余度 |
| `Diagnostics.VarianceComponentIterations` | 方差分量估计迭代次数 |
| `Diagnostics.*Available` | 条件数、完整/站点协方差和残差诊断是否真实可用 |
| `Diagnostics.GlobalTest` | 有自由度且随机模型未由当前数据估计时的双侧卡方整体模型检验；启用方差分量估计或实际发生鲁棒降权时省略 |
| `Warnings` | 条件数过高、鲁棒迭代未收敛或启用自由网内部基准等非致命提示 |

残差符号统一为：

```text
Residual = Observed - Adjusted
Adjusted = Position(To) - Position(From)
```

`AdjustedStation.FormalCovariance` 与 `AdjustedStation.Covariance` 是该站对应的 `3×3` 对角块；跨站相关项请从结果顶层完整矩阵读取。

## JSON 数据交换

`model` 中的输入输出类型都提供 JSON 标签。例如最小输入：

```json
{
  "name": "line",
  "stations": [
    {"id": "A", "fixed": false},
    {"id": "B", "fixed": false}
  ],
  "baselines": [
    {
      "id": "AB",
      "from": "A",
      "to": "B",
      "vector": {"east": 10, "north": 2, "up": 0.5},
      "covariance": {
        "data": [0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0004]
      }
    }
  ],
  "priors": [
    {
      "id": "control-a",
      "station_id": "A",
      "position": {"east": 0, "north": 0, "up": 0},
      "covariance": {
        "data": [0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0004]
      }
    }
  ]
}
```

## 错误判断

公共错误支持 `errors.Is`：

```go
switch {
case errors.Is(err, network.ErrDisconnectedNetwork):
    // 外部基准模式下某个分量没有基准，或自由网模式下存在无观测孤立站
case errors.Is(err, network.ErrInvalidCovariance):
    // 协方差不合法
case errors.Is(err, network.ErrRankDeficient):
    // 数值上仍然秩亏
case errors.Is(err, network.ErrNotConverged):
    // PCG 在指定迭代次数内未达到残差容差
case errors.Is(err, network.ErrInsufficientRedundancy):
    // 某个基线组没有足够的独立冗余度用于估计协方差倍率
case errors.Is(err, network.ErrInvalidProblem):
    // 其他输入校验错误
}
```

字段级校验错误可用 `errors.As` 取得 `*network.ValidationError`，读取 `Field`、`ID` 和 `Message`。

## 实现流程

`network` 包内部按以下顺序工作：

1. 校验站点、基线、坐标先验、有限值和 `3×3` 协方差；
2. 建立无向拓扑，确认每个连通分量都有外部基准，或按配置为无基准分量建立零质心内部基准；
3. 为每个自由站建立 `E/N/U` 三个参数；
4. 按 `-Position(From) + Position(To) = Observed` 组装设计矩阵；
5. 固定站坐标移入观测方程常数项；
6. 组装带完整协方差的控制点坐标先验；
7. 使用完整基线及先验协方差组装稠密或 CSR 稀疏法方程，并采用直接法、PCG 或精确约束零空间投影 PCG 求解；
8. 可选地按 `Baseline.Group` 的目标函数贡献和有效冗余度迭代估计协方差倍率；
9. 可选地按基线 Mahalanobis 距离执行 Huber 迭代重加权；
10. 按协方差模式计算完整逆矩阵、按需局部块或仅整体统计，并生成站点、基线和先验结果。

## 当前边界

- 只处理已经统一到共同 ENU 坐标系的相对向量，不负责经纬高/ECEF/ENU 转换；
- 不负责接收机原始观测、星历、差分、模糊度或 RTK 解算；
- 固定站按精确坐标处理；有不确定度的控制点应通过 `Priors` 提供；
- 支持一条基线内部 E/N/U 相关，不支持不同基线之间的交叉协方差；
- 默认稠密模式适合中小站网；大型站网可选稀疏 PCG。ENU 网络默认使用站点 `3×3` 块 Jacobi，也可选标量 Jacobi；极弱网形、长链或病态站网仍可能需要较多迭代；
- 精确约束 PCG 需要稠密分解约束 Gram 矩阵，适合约束数远少于参数数的场景；若业务创建了大量全局精确约束，其约束层内存和计算量会按约束数增长；
- 方差分量估计当前以完整 ENU 基线组为单位，只估计各组协方差的共同尺度，不分别估计 E/N/U 分量，也不估计控制点先验尺度；
- 自由网当前支持每个无基准连通分量的零质心内部基准，不包含拟稳基准、最小迹基准或求解后的 S 变换；
- 尚未实现速度网和多历元动态模型。

`core/`、`batch/`、`nonlinear/`、`robust/`、`variance/`、`quality/` 和 `sequential/` 保留了底层数值能力，供高级用法扩展；面向站间 ENU 向量网时应优先使用 `model + network` API。

## 验证

```bash
go test ./...
go test -race ./...
go vet ./...
go run ./examples/enu-network
go test ./network -run=^$ -bench=BenchmarkSparsePCGChain500 -benchmem
go test ./network -run=^$ -fuzz=FuzzValidateAndSolveMinimalENU
```

测试覆盖闭合站网、非零固定站、全固定站网、软控制点基准、固定站与控制点混合、相关先验协方差、零质心自由网、无基准连通分量、非法协方差、输入不可变、JSON 往返、整条基线鲁棒降权、分组方差尺度估计、精确约束投影 PCG、块 Jacobi、上下文取消与大型稀疏自由网；同时提供大型链网 benchmark 和最小 ENU 网络 fuzz 入口。

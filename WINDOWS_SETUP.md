# Windows + Visual Studio 配置指南

## 快速开始（3步）

### 第1步：确认你的 Houdini 路径

打开文件资源管理器，找到你的 Houdini 安装目录，通常是：
```
C:\Program Files\Side Effects Software\Houdini 21.0.741
```
（版本号可能不同，找到你实际安装的那个）

### 第2步：安装 Eigen

**推荐方式 A：直接放入项目（最简单，不需要安装任何东西）**

1. 下载 Eigen 3.4.0：https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
2. 解压后，把 `eigen-3.4.0` 文件夹放入项目的 `third_party/` 目录：
   ```
   apic_fluid/
   └── third_party/
       └── eigen-3.4.0/      ← 放这里
           └── Eigen/
               └── Core       ← 确认这个文件存在
   ```
3. CMake 会自动检测到，无需其他配置。

**推荐方式 B：用 vcpkg（如果你已经有 vcpkg）**
```powershell
vcpkg install eigen3:x64-windows
```
VS 和 CMake 会自动找到，无需额外配置。

### 第3步：用 VS 打开项目

1. 打开 Visual Studio 2019 或 2022
2. 菜单：**文件 → 打开 → 文件夹**
3. 选择 `apic_fluid` 文件夹
4. VS 会自动检测 `CMakeLists.txt` 并开始配置

---

## 修改 Houdini 路径

打开 `CMakeSettings.json`（VS 里直接双击），找到这一行：

```json
"value": "C:/Program Files/Side Effects Software/Houdini 21.0.741"
```

**把 `21.0.741` 改成你实际的版本号**（两处，Release 和 Debug 各一处）。

**版本号在哪看？**  
打开 Houdini → Help → About Houdini，版本号如 `21.0.741`。

或者直接在 VS 的 CMake 变量编辑器里改：
- 菜单：**项目 → CMake 设置**
- 找到 `HOUDINI_INSTALL_DIR`，点击修改

---

## 构建步骤

1. VS 打开 folder 后，等待右下角 CMake 配置完成（进度条消失）
2. 菜单：**生成 → 全部生成**（或 Ctrl+Shift+B）
3. 生成的 `SOP_ApicFluid.dll` 会在 `build/x64-Release/` 里

**安装到 Houdini：**
```
菜单：生成 → 安装 apic_fluid
```
这会把 `.dll` 自动复制到 `%USERPROFILE%\Documents\houdini21.0\dso\`

---

## 运行单元测试

VS 里：**测试 → 测试资源管理器**（如果看不到，菜单：视图 → 测试资源管理器）

会看到三个测试：
- `p2g_conservation` — 质量/动量守恒
- `angular_momentum` — APIC vs PIC 角动量
- `pressure_solver` — 压力求解收敛

**注意：单元测试不需要 Houdini，可以独立运行。**

---

## 常见问题

### Q: CMake 提示"找不到 Houdini HDK"
**A:** `CMakeSettings.json` 里的版本号和你实际安装的不一致。  
打开 `C:\Program Files\Side Effects Software\` 看实际文件夹名，复制过来。

### Q: CMake 提示"找不到 Eigen3"
**A:** 按第2步装 Eigen。如果用了 vcpkg，确认 VS 的 vcpkg 集成已开启：  
`vcpkg integrate install`

### Q: 编译报错 `NOMINMAX` 或 `min/max` 相关
**A:** CMakeLists.txt 里已经定义了 `NOMINMAX`，这个错误不应出现。  
如果还是出现，在 `CMakeSettings.json` 里加：
```json
{ "name": "CMAKE_CXX_FLAGS", "value": "/DNOMINMAX" }
```

### Q: 找不到 `SOP_ApicFluid.dll`
**A:** 先确认编译成功（无红色错误）。dll 在：
`build\x64-Release\SOP_ApicFluid.dll`

### Q: Houdini 里看不到 APIC Fluid Solver 节点
**A:** 确认 dll 在 dso 目录里，重启 Houdini。  
dso 目录：`C:\Users\<你的用户名>\Documents\houdini21.0\dso\`

---

## 目录结构

```
apic_fluid/
├── CMakeLists.txt          ← 构建配置（已自动检测路径）
├── CMakeSettings.json      ← VS Open Folder 配置（改这里的版本号）
├── .vs/launch.vs.json      ← 测试启动配置
├── WINDOWS_SETUP.md        ← 本文件
├── README.md               ← 完整文档
├── include/                ← 头文件
├── src/                    ← 核心实现
│   ├── ApicSolver.cpp      ← P2G / G2P / 压力投影
│   ├── ApicGrid.cpp        ← 网格 + B-spline 插值
│   ├── PressureSolver.cpp  ← PCG (Eigen) / Jacobi
│   └── ...
├── tests/                  ← 独立单元测试（无需 Houdini）
├── python/                 ← HDA Python 脚本
└── third_party/            ← 放 eigen-3.4.0 的地方
    └── eigen-3.4.0/        ← 下载解压到这里
```

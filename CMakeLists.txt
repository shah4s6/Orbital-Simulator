cmake_minimum_required(VERSION 3.14)
project(OrbitalSimulator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# ====================================================
#  Fetch project includes
# ====================================================
FetchContent_Declare(imgui GIT_REPOSITORY https://github.com/ocornut/imgui.git GIT_TAG v1.92.7)
FetchContent_Declare(raylib GIT_REPOSITORY https://github.com/raysan5/raylib.git GIT_TAG v5.5)
FetchContent_Declare(rlImGui GIT_REPOSITORY https://github.com/raylib-extras/rlImGui.git GIT_TAG master)

FetchContent_MakeAvailable(imgui raylib rlImGui)

add_executable(OrbitalSimulator src/main.cpp)

target_link_libraries(OrbitalSimulator PRIVATE raylib rlImGui imgui)

target_include_directories(OrbitalSimulator PRIVATE ${imgui_SOURCE_DIR} ${rlImGui_SOURCE_DIR})

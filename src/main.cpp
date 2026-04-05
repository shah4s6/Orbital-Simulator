#include <cmath>
#include "raylib.h"
#include "raymath.h"
#include "../imgui/imgui.h"
#include "../rlImGui/rlImGui.h"

int main()
{
    InitWindow(800, 600, "Orbit Simulator");

    rlImGuiSetup(true);

    Camera3D camera;
    camera.position = (Vector3){10.0f, 10.0f, 10.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 60.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    // Initialize the Earth
    Vector3 earthPos = {0.0f, 0.0f, 0.0f};
    float earthRadius = 1.0f;
    
    // Initialize a satellite
    Vector3 satellitePos = {0.0, 0.0f, 0.0f};
    float satelliteRadius = 0.5f;
    float orbitRadius = 3.5f;
    float orbitSpeed = 1.0f;
    float orbitAngle = 0.0f;

    // UI
    bool show_controls = true;
    bool auto_rotate = true;

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();
        orbitAngle += orbitSpeed * dt;

        if (orbitAngle > 2 * PI) orbitAngle -= 2 * PI;

        satellitePos.x = orbitRadius * cos(orbitAngle);
        satellitePos.y = 0.0f;
        satellitePos.z = orbitRadius * sin(orbitAngle);

        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
        {
            UpdateCamera(&camera, CAMERA_ORBITAL);
        }

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(camera);
        DrawSphere(earthPos, earthRadius, BLUE);
        DrawSphere(satellitePos, satelliteRadius, RED);
        DrawLine3D(earthPos, satellitePos, YELLOW);
        DrawGrid(10, 1.0f);
        EndMode3D();

        rlImGuiBegin();
        ImGui::Begin("Window");
        ImGui::Text("ImGui is working");
        ImGui::Text("Hold the right mouse button to rotate the view");
        ImGui::End();
        rlImGuiEnd();

        EndDrawing();
    }

    CloseWindow();
    return 0;
}

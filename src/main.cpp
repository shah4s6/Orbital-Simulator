#include <cmath>
#include "raylib.h"
#include "raymath.h"
#include "../imgui/imgui.h"
#include "../rlImGui/rlImGui.h"

struct Telemetry
{
    float distance;
    float speed;
    float trueAnomaly;
    float perigee;
    float apogee;
    float semiMajorAxis;
    float eccentricity;
};

// Helper function to draw an elliptical orbit
void DrawOrbitEllipse(Vector3 focus, float semiMajorAxis, float eccentricity, Color color)
{
    if (eccentricity >= 1.0f) return; // e - amound of stretched (0 = circle, <1 = ellipse)

    float c = semiMajorAxis * eccentricity;
    Vector3 center = {focus.x - c, focus.y, focus.z};
    float semiMinorAxis = semiMajorAxis * sqrt(1 - eccentricity * eccentricity);
    
    int segments = 200;
    for (int i = 0; i < segments; i++)
    {
        float angle1 = (float)i / segments * 2 * PI;
        float angle2 = (float)(i + 1) / segments * 2 * PI;
            
        Vector3 p1 = {
            center.x + semiMajorAxis * cos(angle1),
            center.y,
            center.z + semiMinorAxis * sin(angle1)
        };
        Vector3 p2 = {
            center.x + semiMajorAxis * cos(angle2),
            center.y,
            center.z + semiMinorAxis * sin(angle2)
        };
        DrawLine3D(p1, p2, color);
    }
}

float CalculateOrbitalSpeed(float distance, float semiMajorAxis)
{
    // Formula (Keplerian orbit: vis-viva equation): v = sqrt(GM * (2/r - 1/a))
    // v is the relative speed of the two bodies
    // r is the distance between the two bodies' centers
    // a is the length of the semi-major axis
    // G is the gravitational constant
    // M is the mass of the central body
    float velocity = 2.0f / distance - 1.0f / semiMajorAxis;
    return sqrt(fmaxf(0.1f, velocity)) * 3.0f;
}

int main()
{
    InitWindow(800, 600, "Orbit Simulator");
    rlImGuiSetup(true);

    Color DARKERGRAY = GetColor(0x303234FF);

    Camera3D camera;
    camera.position = (Vector3){10.0f, 10.0f, 10.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 60.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    // Initialize the Earth
    Vector3 earthPos = {0.0f, 0.0f, 0.0f};
    float earthRadius = 1.0f;
    
    // Orbit Parameters
    float semiMajorAxis = 5.0f;             // a - half the longest diameter
    float eccentricity = 0.6f;              // e - amound of stretched (0 = circle, <1 = ellipse)
    float orbitSpeed = 1.0f;
    float orbitAngle = 0.0f;

    // Initialize a satellite
    Vector3 satellitePos = {0.0, 0.0f, 0.0f};
    float satelliteRadius = 0.5f;

    Telemetry telemetry = {0};
    bool auto_rotate = true;

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();

        if (auto_rotate)
        {
            orbitAngle += orbitSpeed * dt;
            if (orbitAngle > 2 * PI) orbitAngle -= 2 * PI;
        }

        float c = semiMajorAxis * eccentricity;
        Vector3 ellipseCenter = {earthPos.x - c, earthPos.y, earthPos.z};
        float semiMinorAxis = semiMajorAxis * sqrt(1 - eccentricity * eccentricity);

        // Calculate satellite position
        satellitePos.x = ellipseCenter.x + semiMajorAxis * cos(orbitAngle);
        satellitePos.y = 0.0f;
        satellitePos.z = ellipseCenter.z + semiMinorAxis * sin(orbitAngle);

        // Update telemetry
        telemetry.distance = Vector3Distance(earthPos, satellitePos);
        telemetry.trueAnomaly = orbitAngle * 180.0f / PI;
        telemetry.perigee = semiMajorAxis * (1 - eccentricity);
        telemetry.apogee = semiMajorAxis * (1 + eccentricity);
        telemetry.semiMajorAxis = semiMajorAxis;
        telemetry.eccentricity = eccentricity;
        telemetry.speed = CalculateOrbitalSpeed(telemetry.distance, telemetry.semiMajorAxis);


        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
        {
            UpdateCamera(&camera, CAMERA_ORBITAL);
        }

        BeginDrawing();
        ClearBackground(DARKERGRAY);

        BeginMode3D(camera);
        DrawSphere(earthPos, earthRadius, BLUE);
        DrawSphere(satellitePos, satelliteRadius, RED);
        DrawLine3D(earthPos, satellitePos, YELLOW);
        DrawOrbitEllipse(earthPos, semiMajorAxis, eccentricity, WHITE);
        DrawGrid(10, 1.0f);
        EndMode3D();

        rlImGuiBegin();
        ImGui::Begin("Window");
        ImGui::Text("Hold the right mouse button to rotate the view");
        ImGui::TextColored(ImVec4 (0.0f, 1.0f, 0.0f, 1.0f), "Perigee");
        ImGui::TextColored(ImVec4 (1.0f, 1.0f, 0.0f, 1.0f), "Apogee");
        ImGui::SliderFloat("Semi-Major Axis (a)",&semiMajorAxis, 2.0f, 8.0f);
        ImGui::SliderFloat("Eccentricity (e)",&eccentricity, 0.0f, 0.99f);
        ImGui::SliderFloat("Orbit Speed",&orbitSpeed, 0.0f, 2.0f);
        ImGui::Checkbox("Automatic Rotation", &auto_rotate);
        ImGui::Separator();
        ImGui::Text("Distance: %.2f", telemetry.distance);
        ImGui::Text("Speed: %.2f", telemetry.speed);
        ImGui::Text("True Anomaly: %.1f", telemetry.trueAnomaly);
        ImGui::Text("Perigee: %.2f | Apogee: %.2f", telemetry.perigee, telemetry.apogee);
        ImGui::Separator();

        if (ImGui::Button("Reset"))
        {
            semiMajorAxis = 5.0f;
            eccentricity = 0.6f;
            orbitSpeed = 1.0f;
            orbitAngle = 0.0f;
            auto_rotate = true;
        }

        ImGui::End();
        rlImGuiEnd();

        EndDrawing();
    }

    CloseWindow();
    return 0;
}

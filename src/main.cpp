#include <cmath>
#include "raylib.h"
#include "raymath.h"
#include "../imgui/imgui.h"
#include "../rlImGui/rlImGui.h"
#include "eigen3/Eigen/Dense"

struct Telemetry
{
    float distance;
    float velocity;
    float trueAnomaly;
    float perigee;
    float apogee;
    float semiMajorAxis;
    float eccentricity;
};

struct Satellite
{
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double mass;
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

// RK4 ODE: dy/dt = f(t, y)
// Runge-Kutta 4th Order
// void RK4_Method(double(*f)(double t, double y), double y0, double t0, double h)
// {
//     double t = t0;
//     double y = y0;

//     // Compute the four slopes
//     double k1 = f(t, y);
//     double k2 = f(t + h/2.0, y + (k1/2.0) * h);
//     double k3 = f(t + h/2.0, y + (k2/2.0) * h);
//     double k4 = f(t + h, y + k3 * h);

//     // Update y using weighted average of slopes
//     y = y + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4) * h;

//     return y;
// }

// Satellite has 6 state variables (position and velocity in the x, y, z directions)

struct StateVariables
{
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
};

// Using ODE45 like in MATLAB
StateVariables calculateDerivatives(const StateVariables& state, double t, double GM)
{
    StateVariables derivative;
    // r = position, dr/dt = velocity
    derivative.position = state.velocity;

    // dv/dt = acceleration
    Eigen::Vector3d r = state.position;
    double r_magnitude = r.norm();

    if (r_magnitude > 0.01)
    {
        // F = GM m / r^2 = ma -> a = -GM*r_hat/r^3
        derivative.velocity = -GM * r / (r_magnitude * r_magnitude * r_magnitude);
    }
    else
    {
        derivative.velocity = Eigen::Vector3d(0.0, 0.0, 0.0);
    }

    return derivative;
}

void rk4Method(Satellite& sat, double dt, double GM)
{
    // Initialize the time
    double t = 0.0;

    // Current state
    StateVariables current;
    current.position = sat.position;
    current.velocity = sat.velocity;

    // k1: Derivative at the start of the interval
    StateVariables k1 = calculateDerivatives(current, t, GM);

    // k2: Derivative at the midpoint
    StateVariables state_k2;
    state_k2.position = current.position + k1.position * (dt / 2.0);
    state_k2.velocity = current.velocity + k1.velocity * (dt / 2.0);
    StateVariables k2 = calculateDerivatives(state_k2, t + dt / 2.0, GM);

    // k3: Derivative at another midpoint
    StateVariables state_k3;
    state_k3.position = current.position + k2.position * (dt / 2.0);
    state_k3.velocity = current.velocity + k2.velocity * (dt / 2.0);
    StateVariables k3 = calculateDerivatives(state_k3, t + dt / 2.0, GM);

    // k4: Derivative at the end of the interval
    StateVariables state_k4;
    state_k4.position = current.position + k3.position * dt;
    state_k4.velocity = current.velocity + k3.velocity * dt;
    StateVariables k4 = calculateDerivatives(state_k4, t + dt, GM);

    // y_new = y_current + (k1 + 2k2 + 2k3 + k4) * dt/6
    sat.position = current.position + (k1.position + 2.0*k2.position + 2.0*k3.position + k4.position) * (dt / 6.0);
    sat.velocity = current.velocity + (k1.velocity + 2.0*k2.velocity + 2.0*k3.velocity + k4.velocity) * (dt / 6.0);

}

int main()
{
    InitWindow(1080, 720, "Orbit Simulator");
    rlImGuiSetup(true);

    const double GM = 15.0;

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

    // Initialize a satellite
    Satellite sat;
    double r_initial = semiMajorAxis * (1 - eccentricity);
    sat.position = Eigen::Vector3d(r_initial, 0.0, 0.0);
    double v_perigee = sqrt(GM * (1 + eccentricity) / (semiMajorAxis * (1 - eccentricity)));
    sat.velocity = Eigen::Vector3d(0.0, 0.0, v_perigee);
    sat.mass = 1.0;
    float satelliteRadius = 0.5f;

    Telemetry telemetry = {0};
    bool auto_rotate = true;

    bool use2p5D = true;

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();

        if (auto_rotate)
        {
            Eigen::Vector3d r = sat.position;
            double r_magnitude = r.norm();

            Eigen::Vector3d acceleration = -GM * r / (r_magnitude * r_magnitude * r_magnitude);

            rk4Method(sat, dt, GM);

            if (use2p5D)
            {
                sat.position.y() = 0.0;
                sat.velocity.y() = 0.0;
            }
        }

        if (use2p5D)
        {
            // Fixed camera for debugging
            camera.position = (Vector3){0.0f, 15.0f, 10.0f};
            camera.target = (Vector3){0.0f, 0.0f, 0.0f};
            camera.up = (Vector3){0.0f, 1.0f, 0.0f};
        }
        else if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
        {
            UpdateCamera(&camera, CAMERA_FIRST_PERSON);
        }
        
        // Update telemetry
        Eigen::Vector3d r = sat.position;
        Eigen::Vector3d v = sat.velocity;
        double r_magnitude = r.norm();
        double v_magnitude = v.norm();

        telemetry.distance = r_magnitude;
        telemetry.velocity = v_magnitude;

        double epsilon = v_magnitude * v_magnitude / 2.0 - GM / r_magnitude; // epsilon - specific orbital energy 
        telemetry.semiMajorAxis = -GM / (2.0 * epsilon);

        // Eccentricity vector
        Eigen::Vector3d e_vector = ((v_magnitude * v_magnitude - GM / r_magnitude) * r - (r.dot(v)) * v) / GM;
        telemetry.eccentricity = e_vector.norm();

        // Perigee and Apogee
        telemetry.perigee = telemetry.semiMajorAxis * (1 - telemetry.eccentricity);
        telemetry.apogee = telemetry.semiMajorAxis * (1 + telemetry.eccentricity);

        telemetry.trueAnomaly = acos(e_vector.dot(r) / (e_vector.norm() * r_magnitude)) * 180.0 / PI;
        if (r.dot(v) < 0) telemetry.trueAnomaly = 360.0 - telemetry.trueAnomaly;

        // Calculate perigeee and apogee marker positions
        float perigeeDistance = (float)telemetry.perigee;
        float apogeeDistance = (float)telemetry.apogee;

        Vector3 perigeePos = {earthPos.x + perigeeDistance, earthPos.y, earthPos.z};
        Vector3 apogeePos = {earthPos.x - apogeeDistance, earthPos.y, earthPos.z};

        // Calculate satellite position
        Vector3 satellitePos = {(float)sat.position.x(), (float)sat.position.y(), (float)sat.position.z()};

        BeginDrawing();
        ClearBackground(DARKERGRAY);

        BeginMode3D(camera);
        DrawSphere(earthPos, earthRadius, BLUE);
        DrawSphere(satellitePos, satelliteRadius, RED);
        DrawLine3D(earthPos, satellitePos, YELLOW);
        DrawOrbitEllipse(earthPos, (float)telemetry.semiMajorAxis, (float)telemetry.eccentricity, WHITE);
        DrawSphere(perigeePos, 0.1f, GREEN);
        DrawSphere(apogeePos, 0.1f, YELLOW);
        DrawGrid(10, 1.0f);
        EndMode3D();

        rlImGuiBegin();
        ImGui::Begin("Window");
        ImGui::Text("Hold the right mouse button to rotate the view");

        // ImGui::SliderFloat("Semi-Major Axis (a)",&semiMajorAxis, 2.0f, 8.0f);    TODO: Reimplement this slider correctly
        // ImGui::SliderFloat("Eccentricity (e)",&eccentricity, 0.0f, 0.99f);       TODO: Reimplement this slider correctly
        // ImGui::SliderFloat("Speed", &telemetry.velocity, 0.0f, 5.0f);            TODO:Reimplement this slider correctly
        ImGui::Checkbox("Automatic Rotation", &auto_rotate);
        ImGui::Separator();
        ImGui::Text("Distance: %.2f", telemetry.distance);
        ImGui::Text("True Anomaly: %.1f", telemetry.trueAnomaly);
        ImGui::TextColored(ImVec4 (0.0f, 1.0f, 0.0f, 1.0f), "Perigee: %.2f", telemetry.perigee);
        ImGui::TextColored(ImVec4 (1.0f, 1.0f, 0.0f, 1.0f), "Apogee: %.2f", telemetry.apogee);
        ImGui::Separator();

        ImGui::Checkbox("Debug Mode (Top Down)", &use2p5D);

        if (ImGui::Button("Reset"))
        {
            semiMajorAxis = 5.0f;
            eccentricity = 0.6f;
            auto_rotate = true;
            
            // Reset satellite state
            double r_initial = semiMajorAxis * (1 - eccentricity);
            sat.position = Eigen::Vector3d(r_initial, 0.0, 0.0);
            double v_perigee = sqrt(GM * (1 + eccentricity) / (semiMajorAxis * (1 - eccentricity)));
            sat.velocity = Eigen::Vector3d(0.0, 0.0, v_perigee);
        }

        ImGui::End();
        rlImGuiEnd();

        EndDrawing();
    }

    CloseWindow();
    return 0;
}

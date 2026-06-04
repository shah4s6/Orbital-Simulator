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
    // TODO: Add ground track
};

struct Satellite
{
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double mass;
};

// Helper function to draw an elliptical orbit
void DrawOrbitEllipse(Vector3 focus, float semiMajorAxis, float eccentricity, float inclination, Color color)
{
    if (eccentricity >= 1.0f) return; // e - amound of stretched (0 = circle, <1 = ellipse)

    float c = semiMajorAxis * eccentricity;
    Vector3 center = {focus.x - c, focus.y, focus.z};
    float semiMinorAxis = semiMajorAxis * sqrt(1 - eccentricity * eccentricity);

    int segments = 200;
    float inclinationRadian = inclination * PI / 180.0f;
    float cosInclination = cos(inclinationRadian);
    float sinInclination = sin(inclinationRadian);
    for (int i = 0; i < segments; i++)
    {
        float angle1 = (float)i / segments * 2 * PI;
        float angle2 = (float)(i + 1) / segments * 2 * PI;

        Vector3 p1 = {
            center.x + semiMajorAxis * cos(angle1),
            center.y + semiMinorAxis * sin(angle1) * sinInclination,
            center.z + semiMinorAxis * sin(angle1) * cosInclination
        };
        Vector3 p2 = {
            center.x + semiMajorAxis * cos(angle2),
            center.y + semiMinorAxis * sin(angle2) * sinInclination,
            center.z + semiMinorAxis * sin(angle2) * cosInclination
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

float CalculateInclination(const Eigen::Vector3d& position, const Eigen::Vector3d & velocity)
{
    Eigen::Vector3d angularMomentum = position.cross(velocity);
    double incline = acos(angularMomentum.z() / angularMomentum.norm()) * 180.0 / PI;
    return incline;
}

float CalculateRequiredDeltaV(double currentSemiMajorAxis, double targetSemiMajorAxis, double GM, double currentVel)
{
    // Hohmann transfer approximation
    double currentOrbitSpeed = sqrt(GM / currentSemiMajorAxis);
    double targetOrbitSpeed = sqrt(GM / targetSemiMajorAxis);
    return targetOrbitSpeed - currentOrbitSpeed;
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
    double v_mag = sat.velocity.norm();
    if (v_mag > 20.0) 
    {
        sat.velocity = sat.velocity.normalized() * 20.0;
    }

    double r_mag = sat.position.norm();
    if (r_mag > 30.0) 
    {
        // If too far, scale velocity to bring it back
        sat.velocity *= 0.95;
    }

    // Initialize the time
    double t = 0.0;

    if (dt > 0.1) dt = 0.1;

    if (sat.position.norm() < 0.5) {
        // Too close to Earth - don't integrate further
        return;
    }

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

void ApplyBurn(Satellite& sat, double deltaV, const Eigen::Vector3d& direction, double GM, double dt)
{
    sat.velocity += direction * deltaV;

    Eigen::Vector3d r = sat.position;
    Eigen::Vector3d v = sat.velocity;
    double r_magnitude = r.norm();
    double v_magnitude = v.norm();

    double epsilon = v_magnitude * v_magnitude / 2.0 - GM / r_magnitude;
    double semiMajorAxis = -GM / (2.0 * epsilon);

    Eigen::Vector3d e_vector = ((v_magnitude * v_magnitude - GM / r_magnitude) * r - (r.dot(v)) * v) / GM;
    double eccentricity = e_vector.norm();
}

struct ManueverState
{
    bool isActive;
    float duration;
    float progress;
    Eigen::Vector3d deltaVPerSecond;
    std::string name;
};

int main()
{
    InitWindow(1400, 1040, "Orbit Simulator");
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
    float satelliteRadius = 0.12f;

    Telemetry telemetry = {0};

    // Trail for the satellite path
    std::vector<Vector3> trail;
    float trailTimer = 0.0f;

    // UI
    bool auto_rotate = true;
    bool use2p5D = true;
    static float semiMajorAxisUI = semiMajorAxis;
    static float eccentricityUI = eccentricity;
    bool orbitChanged = false;
    bool isManeuvering = false;
    float maneuverProgress = 0.0f;
    bool perigeeBurn = false;
    bool apogeeBurn = false;
    static bool needsReset = false;
    static float lastSemiMajorAxis = semiMajorAxisUI;
    static float lastEccentricity = eccentricityUI;

    // Target orbit visualization
    float targetSemiMajorAxis = semiMajorAxis;
    float targetEccentricity = eccentricity;
    float targetDisplayInclination = 0.0f;
    bool showTargetOrbit = false;

    // Inclination variable to enable 3D orbits
    float currentInclination = 0.0f;           // Current inclination (degrees)
    float targetInclination = 0.0f;            // Target inclination obtained by the slider (degrees)
    static float timeWarp = 1.0f; // Initialized at normal speed 1x
    float targetOrbitTimer = 0.0f;

    ManueverState activeManeuver = {false, 0.0f, 0.0f, Eigen::Vector3d::Zero(), ""};
    static float burnAmount = 0.5f;
    static float maneuverDuration = 1.0f;

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();
        if (dt > 0.033f) dt = 0.033f; // Capping the dt to prevent skipping of the trajectory

        if (IsWindowMinimized())
        {
            dt = 0.0f;
        }

        if (auto_rotate)
        {
            // FIRST: Apply any active maneuver delta-V for this frame
            if (activeManeuver.isActive)
            {
                activeManeuver.progress += dt * timeWarp;
                float t = activeManeuver.progress / activeManeuver.duration;
                
                if (t >= 1.0f)
                {
                    // Apply remaining delta-V
                    float remainingFraction = 1.0f - ((activeManeuver.progress - dt * timeWarp) / activeManeuver.duration);
                    if (remainingFraction > 0.01f) {
                        sat.velocity += activeManeuver.deltaVPerSecond * remainingFraction;
                    }
                    activeManeuver.isActive = false;
                    TraceLog(LOG_INFO, "Maneuver completed");
                }
                else
                {
                    // Apply delta-V for this frame
                    sat.velocity += activeManeuver.deltaVPerSecond * (dt * timeWarp / activeManeuver.duration);
                }
            }
            
            // THEN: Integrate physics with the updated velocity
            float dtNew = dt * timeWarp;
            rk4Method(sat, dtNew, GM);
        }

        // Remove the old maneuver application code that was here
        
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

        Vector3 perigeePos = {0,0,0};
        Vector3 apogeePos = {0,0,0};

        float maxDisplayDistance = 20.0f;  // Don't draw beyond this range
        float clampedPerigee = fminf(perigeeDistance, maxDisplayDistance);
        float clampedApogee = fminf(apogeeDistance, maxDisplayDistance);
        
        // if (telemetry.eccentricity > 0.01f)  // Non-circular orbit
        // {
        //     // Use eccentricity vector direction
        //     Vector3 perigeeDir = {(float)e_vector.x(), (float)e_vector.y(), (float)e_vector.z()};
        //     float len = sqrt(perigeeDir.x * perigeeDir.x + perigeeDir.y * perigeeDir.y + perigeeDir.z * perigeeDir.z);
        //     if (len > 0.01f)
        //     {
        //         Vector3 norm = {perigeeDir.x / len, perigeeDir.y / len, perigeeDir.z / len};
        //         perigeePos = {earthPos.x + perigeeDistance * norm.x, earthPos.y + perigeeDistance * norm.y, earthPos.z + perigeeDistance * norm.z};
        //         apogeePos = {earthPos.x - apogeeDistance * norm.x, earthPos.y - apogeeDistance * norm.y, earthPos.z - apogeeDistance * norm.z};
        //     }
        // }
        // else  // Circular orbit - use radial direction
        // {
        //     Vector3 radialDir = {(float)r.x(), (float)r.y(), (float)r.z()};
        //     float len = sqrt(radialDir.x * radialDir.x + radialDir.y * radialDir.y + radialDir.z * radialDir.z);
        //     if (len > 0.01f)
        //     {
        //         Vector3 norm = {radialDir.x / len, radialDir.y / len, radialDir.z / len};
        //         perigeePos = {earthPos.x + perigeeDistance * norm.x, earthPos.y + perigeeDistance * norm.y, earthPos.z + perigeeDistance * norm.z};
        //         apogeePos = {earthPos.x - apogeeDistance * norm.x, earthPos.y - apogeeDistance * norm.y, earthPos.z - apogeeDistance * norm.z};
        //     }
        // }

        if (telemetry.eccentricity >= 0.0f && telemetry.eccentricity < 0.99f && 
            telemetry.semiMajorAxis > 0.1f && telemetry.semiMajorAxis < 100.0f)
        {
            // Get the eccentricity vector direction (points to perigee)
            if (e_vector.norm() > 0.01f)
            {
                Vector3 perigeeDir = {(float)e_vector.x(), (float)e_vector.y(), (float)e_vector.z()};
                float len = Vector3Length(perigeeDir);
                if (len > 0.01f)
                {
                    Vector3 norm = {perigeeDir.x / len, perigeeDir.y / len, perigeeDir.z / len};
                    
                    // Perigee is in direction of eccentricity vector
                    perigeePos = {earthPos.x + clampedPerigee * norm.x, 
                                earthPos.y + clampedPerigee * norm.y, 
                                earthPos.z + clampedPerigee * norm.z};
                    // Apogee is opposite direction
                    apogeePos = {earthPos.x - clampedApogee * norm.x, 
                                earthPos.y - clampedApogee * norm.y, 
                                earthPos.z - clampedApogee * norm.z};
                }
            }
        }

        // Calculate satellite position
        Vector3 satellitePos = {(float)sat.position.x(), (float)sat.position.y(), (float)sat.position.z()};

        // Update trail
        trailTimer += dt;
        float trailInterval = 0.05 / timeWarp;
        if (trailTimer >= trailInterval)
        {
            trailTimer = 0.0f;
            trail.push_back(satellitePos);

            while (trail.size() > 800)
            {
                trail.erase(trail.begin());
            }
        }

        BeginDrawing();
        ClearBackground(DARKERGRAY);

        BeginMode3D(camera);
        DrawSphere(earthPos, earthRadius, DARKBLUE);
        DrawSphereWires(earthPos, earthRadius, 16, 16, LIGHTGRAY);
        DrawSphere(satellitePos, satelliteRadius, RED);
        
        // Draw a velocity vector
        if (v_magnitude > 0.01f)
        {
            Vector3 start = satellitePos;
            Vector3 direction = {(float)v.x(), (float)v.y(), (float)v.z()};

            float invMagnitude = 1.0f / v_magnitude;
            direction.x *= invMagnitude;
            direction.y *= invMagnitude;
            direction.z *= invMagnitude;

            float arrowLength = fminf(v_magnitude * 0.5f, 2.0f);
            Vector3 end = {start.x + direction.x * arrowLength, start.y + direction.y * arrowLength, start.z + direction.z * arrowLength};

            DrawLine3D(start, end, ORANGE);
        }

        DrawLine3D(earthPos, satellitePos, YELLOW);
        DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){10.0f, 0.0f, 0.0f}, RED);
        DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 10.0f, 0.0f}, GREEN);
        DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 0.0f, 10.0f}, BLUE);

        if (showTargetOrbit)
        {
            // Draw ghost target orbit when sliders are adjusted
            Color targetColor = GREEN;

            targetOrbitTimer += dt;
            if (targetOrbitTimer > 3.0f)
            {
                showTargetOrbit = false;
                targetOrbitTimer = 0.0f;
            }
        }
        else
        {
            targetOrbitTimer = 0.0f;
        }

        if (trail.size() > 1)
        {
            for (size_t i = 0; i < trail.size() - 1; i++)
            {
                Color trailColor = MAGENTA;
                DrawLine3D(trail[i], trail[i + 1], trailColor);
            }
        }

        DrawSphere(perigeePos, 0.1f, GREEN);
        DrawSphere(apogeePos, 0.1f, YELLOW);
        DrawGrid(10, 1.0f);
        EndMode3D();

        rlImGuiBegin();
        // Telemetry Panel (Left)
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(260, 300), ImGuiCond_Always);
        ImGui::Begin("Telemetry", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        ImGui::Text("Orbital Data");
        ImGui::Separator();
        
        ImGui::Text("Distance: %.2f", telemetry.distance);
        ImGui::Text("Speed: %.2f", telemetry.velocity);
        ImGui::Text("True Anomaly: %.2f", telemetry.trueAnomaly);
        ImGui::Separator();

        ImGui::Text("Orbital Elements");
        ImGui::Text("Semi-Major Axis(a): %.2f", telemetry.semiMajorAxis);
        ImGui::Text("Eccentricity: %.2f", telemetry.eccentricity);
        ImGui::Separator();

        ImGui::TextColored(ImVec4 (0.0f, 1.0f, 0.0f, 1.0f), "Perigee: %.2f", telemetry.perigee);
        ImGui::TextColored(ImVec4 (1.0f, 1.0f, 0.0f, 1.0f), "Apogee: %.2f", telemetry.apogee);
        ImGui::Separator();
        
        ImGui::TextColored(ImVec4(0.502, 0.502, 0.502, 1.0), "1 unit = 6371 km");

        ImGui::End();

        // Orbit Controls Panel (Right)
        ImGui::SetNextWindowPos(ImVec2(GetScreenWidth() - 270, 10), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(260, 320), ImGuiCond_Always);
        ImGui::Begin("Orbit Controls", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        ImGui::Text("Parameters");
        ImGui::Separator();

        ImGui::Text("Semi-Major Axis (a)");
        if (ImGui::SliderFloat("##semiMajorAxis", &semiMajorAxisUI, 2.0f, 8.0f))
        {
            orbitChanged = true;
            targetSemiMajorAxis = semiMajorAxisUI;
            showTargetOrbit = true;
        }
        ImGui::Text("Eccentricity (e)");
        if (ImGui::SliderFloat("##eccentricity", &eccentricityUI, 0.0f, 0.99f))
        {
            orbitChanged = true;
            targetEccentricity = eccentricityUI;
            showTargetOrbit = true;
        }
        ImGui::Text("Inclination (i)");
        if (ImGui::SliderFloat("##inclination", &targetInclination, 0, 180, "%0.0f"))
        {
            targetDisplayInclination = targetInclination;
            showTargetOrbit = true;
        }

        // Display orbit type based on the inclination
        if (targetInclination < 90.0)
        {
            ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f), "PROGRADE (i < 90°)");
        }
        else if (targetInclination > 90.0f)
        {
            ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.5f, 1.0f), "RETROGRADE (i > 90°)");
        }
        else
        {
            ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1.0f), "POLAR (i = 90°)");
        }
        ImGui::Separator();

        if (showTargetOrbit)
        {
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Target Orbit:");
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Semi-Major Axis %.2f", targetSemiMajorAxis);
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Eccentricity: %.2f", targetEccentricity);
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Inclination: %.0f", targetDisplayInclination);
        }

        if (ImGui::Button("Reset Orbit"))
        {
            semiMajorAxis = 5.0f;
            eccentricity = 0.6f;
            semiMajorAxisUI = 5.0f;
            eccentricityUI = 0.6f;
            currentInclination = 0.0f;
            targetInclination = 0.0f;
            targetDisplayInclination = 0.0f;
            targetSemiMajorAxis = 5.0f;
            targetEccentricity = 0.6f;
            auto_rotate = true;
            showTargetOrbit = false;
            maneuverProgress = 0.0f;
            activeManeuver.isActive = false;
            isManeuvering = false;
            maneuverProgress = 0.0f;

            trail.clear();

            for (int i = 0; i < 5; i++) 
            {
                rk4Method(sat, 0.01, GM);
            }
            
            // Reset satellite state
            double r_initial = semiMajorAxis * (1 - eccentricity);
            sat.position = Eigen::Vector3d(r_initial, 0.0, 0.0);
            double v_perigee = sqrt(GM * (1 + eccentricity) / (semiMajorAxis * (1 - eccentricity)));
            sat.velocity = Eigen::Vector3d(0.0, 0.0, v_perigee);
        }
        ImGui::End();

        // Maneuver Controls Panel (Bottom Right)
        // ImGui::SetNextWindowPos(ImVec2(GetScreenWidth() - 270, GetScreenHeight() - 255), ImGuiCond_Always);
        // ImGui::SetNextWindowSize(ImVec2(260, 245), ImGuiCond_Always);
        // ImGui::Begin("Maneuvers", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        // ImGui::Text("Manual Burns");
        // ImGui::Separator();

        // if (ImGui::Button("Burn at Perigee"))
        // {
        //     perigeeBurn = true;
        // }
        // if (ImGui::Button("Burn at apogee"))
        // {
        //     apogeeBurn = true;
        // }

        // if (perigeeBurn || apogeeBurn)
        // {
        //     if (ImGui::Button("Cancel"))
        //     {
        //         perigeeBurn = false;
        //         apogeeBurn = false;
        //     }
        // }

        // static float burnAmount = 0.5f;
        // ImGui::SliderFloat("Delta-V", &burnAmount, 0.1f, 2.0f, "%.1f m/s");

        // ImGui::Text("Prograde/Retrograde");
        // if (ImGui::Button("+ Prograde"))
        // {
        //     Eigen::Vector3d v_direction = sat.velocity.normalized();
        //     sat.velocity += v_direction * burnAmount;
        // }
        // if (ImGui::Button("- Retrograde"))
        // {
        //     Eigen::Vector3d v_direction = sat.velocity.normalized();
        //     sat.velocity -= v_direction * burnAmount;
        // }

        // ImGui::End();

        // Maneuver Controls Panel (Simplified)
        ImGui::SetNextWindowPos(ImVec2(GetScreenWidth() - 270, GetScreenHeight() - 250), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(260, 240), ImGuiCond_Always);
        ImGui::Begin("Maneuvers", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        ImGui::Text("Manual Burns");
        ImGui::Separator();

        ImGui::SliderFloat("Delta-V", &burnAmount, 0.1f, 3.0f, "%.1f units/s");
        ImGui::SliderFloat("Burn Duration", &maneuverDuration, 0.5f, 3.0f, "%.1f sec");

        ImGui::Separator();
        ImGui::Text("Prograde / Retrograde");

        // Prograde burn
        if (ImGui::Button("+ Prograde"))
        {
            if (!activeManeuver.isActive && sat.velocity.norm() > 0.01)
            {
                Eigen::Vector3d v_direction = sat.velocity.normalized();
                activeManeuver.deltaVPerSecond = v_direction * burnAmount;
                activeManeuver.duration = maneuverDuration;
                activeManeuver.progress = 0.0f;
                activeManeuver.isActive = true;
                activeManeuver.name = "Prograde";
                TraceLog(LOG_INFO, "Starting Prograde burn");
            }
        }
        ImGui::SameLine();

        // Retrograde burn
        if (ImGui::Button("- Retrograde"))
        {
            if (!activeManeuver.isActive && sat.velocity.norm() > 0.01)
            {
                Eigen::Vector3d v_direction = sat.velocity.normalized();
                activeManeuver.deltaVPerSecond = -v_direction * burnAmount;
                activeManeuver.duration = maneuverDuration;
                activeManeuver.progress = 0.0f;
                activeManeuver.isActive = true;
                activeManeuver.name = "Retrograde";
                TraceLog(LOG_INFO, "Starting Retrograde burn");
            }
        }

        ImGui::Separator();
        ImGui::Text("Radial Burns (Small Amounts!)");

        // Radial out (away from Earth) - MUCH smaller default
        if (ImGui::Button("Radial Out"))
        {
            if (!activeManeuver.isActive)
            {
                // Get the direction perpendicular to velocity (radial direction)
                Eigen::Vector3d radial_dir = sat.position.normalized();
                
                // Use a much smaller delta-v for radial burns
                float radialBurnAmount = burnAmount * 0.15f;  // 15% of normal burn
                
                activeManeuver.deltaVPerSecond = radial_dir * radialBurnAmount;
                activeManeuver.duration = maneuverDuration * 0.5f;  // Shorter duration
                activeManeuver.progress = 0.0f;
                activeManeuver.isActive = true;
                activeManeuver.name = "Radial Out";
                TraceLog(LOG_INFO, "Starting Radial Out burn (%.2f units/s)", radialBurnAmount);
            }
        }
        ImGui::SameLine();

        // Radial in (toward Earth)
        if (ImGui::Button("Radial In"))
        {
            if (!activeManeuver.isActive)
            {
                Eigen::Vector3d radial_dir = sat.position.normalized();
                
                // Use a much smaller delta-v for radial burns
                float radialBurnAmount = burnAmount * 0.15f;  // 15% of normal burn
                
                activeManeuver.deltaVPerSecond = -radial_dir * radialBurnAmount;
                activeManeuver.duration = maneuverDuration * 0.5f;
                activeManeuver.progress = 0.0f;
                activeManeuver.isActive = true;
                activeManeuver.name = "Radial In";
                TraceLog(LOG_INFO, "Starting Radial In burn (%.2f units/s)", radialBurnAmount);
            }
        }

        // Cancel button for active maneuver
        if (activeManeuver.isActive)
        {
            ImGui::Separator();
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "Active: %s", activeManeuver.name.c_str());
            ImGui::ProgressBar(activeManeuver.progress / activeManeuver.duration);
            
            if (ImGui::Button("Cancel"))
            {
                activeManeuver.isActive = false;
                TraceLog(LOG_INFO, "Maneuver cancelled");
            }
        }

        // Help text
        ImGui::Separator();
        ImGui::TextColored(ImVec4(0.6f, 0.6f, 0.6f, 1.0f), "Prograde: Raises orbit");
        ImGui::TextColored(ImVec4(0.6f, 0.6f, 0.6f, 1.0f), "Retrograde: Lowers orbit");
        ImGui::TextColored(ImVec4(0.6f, 0.6f, 0.6f, 1.0f), "Radial: Changes eccentricity");

        ImGui::End();


        // Time Controls Panel (Bottom)
        ImGui::SetNextWindowPos(ImVec2(GetScreenWidth()/2 - 200, GetScreenHeight() - 150), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(400, 140), ImGuiCond_Always);
        ImGui::Begin("Time Controls", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        ImGui::Separator();

        ImGui::SliderFloat("Time Warp", &timeWarp, 0.1f, 20.0f, "%.1fx");
        ImGui::SameLine();
        ImGui::Text("Simulation Speed: %.1fx", timeWarp);

        Eigen::Vector3d angularMomentum = sat.position.cross(sat.velocity);
        currentInclination = acos(angularMomentum.z() / angularMomentum.norm()) * 180.0 / PI;

        // if (orbitChanged && auto_rotate)
        // {
        //     semiMajorAxis = semiMajorAxisUI;
        //     eccentricity = eccentricityUI;

        //     double currentTrueAnomalyRadian = telemetry.trueAnomaly * PI / 180.0;
        //     double r_new = semiMajorAxis * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(currentTrueAnomalyRadian);
        //     double v_new = sqrt(GM * (2.0 / r_new - 1/ semiMajorAxis));

        //     Eigen::Vector3d e_vector = ((v_magnitude * v_magnitude - GM / r_magnitude) * r(r.dot(v)) * v) / GM;
        //     Eigen::Vector3d perigeeDirection = e_vector.normalized();
        //     Eigen::Vector3d normalDirection = r.cross(v).normalized();
        //     Eigen::Vector3d perpendicularDirection = normalDirection.cross(perigeeDirection);

        //     double cosTrueAnomaly = cos(currentTrueAnomalyRadian);
        //     double sinTrueAnomaly = sin(currentTrueAnomalyRadian);

        //     sat.position = perigeeDirection * r_new * cosTrueAnomaly + perpendicularDirection * r_new * sinTrueAnomaly;

        //     double h = sqrt(GM * semiMajorAxis * (1- eccentricity * eccentricity));
        // }

        // After the burn application, add:
        if (activeManeuver.isActive == false) {
            // Check for unrealistic velocity
            if (sat.velocity.norm() > 15.0) {
                TraceLog(LOG_WARNING, "Velocity too high (%.2f), limiting", sat.velocity.norm());
                sat.velocity = sat.velocity.normalized() * 15.0;
            }
        }

        if (orbitChanged && auto_rotate)
        {
            semiMajorAxis = semiMajorAxisUI;
            eccentricity = eccentricityUI;

            // Recalculate position at current true anomaly
            double currentTrueAnomalyRadian = telemetry.trueAnomaly * PI / 180.0;
            double r_new = semiMajorAxis * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(currentTrueAnomalyRadian));
            double v_new = sqrt(GM * (2.0 / r_new - 1.0 / semiMajorAxis));

            // Get orbital plane vectors
            Eigen::Vector3d e_vector_calc = ((v_magnitude * v_magnitude - GM / r_magnitude) * r - (r.dot(v)) * v) / GM;
            Eigen::Vector3d perigeeDirection = e_vector_calc.normalized();
            Eigen::Vector3d normalDirection = r.cross(v).normalized();
            Eigen::Vector3d perpendicularDirection = normalDirection.cross(perigeeDirection);

            double cosTrueAnomaly = cos(currentTrueAnomalyRadian);
            double sinTrueAnomaly = sin(currentTrueAnomalyRadian);

            // Update position
            sat.position = perigeeDirection * r_new * cosTrueAnomaly + perpendicularDirection * r_new * sinTrueAnomaly;
            
            // Update velocity magnitude and direction
            double h = sqrt(GM * semiMajorAxis * (1 - eccentricity * eccentricity));
            double velocityRadial = (GM / h) * eccentricity * sin(currentTrueAnomalyRadian);
            double velocityTransverse = (GM / h) * (1 + eccentricity * cos(currentTrueAnomalyRadian));
            
            sat.velocity = perigeeDirection * velocityRadial + perpendicularDirection * velocityTransverse;
        }

        if (orbitChanged && auto_rotate && !isManeuvering)
        {
            // Start manuevering to change the orbit
            isManeuvering = true;
            maneuverProgress = 0.0f;
            orbitChanged = false;
        }

        if (isManeuvering && auto_rotate)
        {
            maneuverProgress += dt * timeWarp;

            float t = fminf(maneuverProgress / 3.0f, 1.0f);

            double targetSemiMajorAxis = semiMajorAxisUI;
            double currentSemiMajorAxis = telemetry.semiMajorAxis;

            if (fabs(currentSemiMajorAxis - targetSemiMajorAxis) > 0.05f)
            {
                Eigen::Vector3d v_direction = sat.velocity.normalized();
                double deltaV = (targetSemiMajorAxis - currentSemiMajorAxis) * 0.05 * dt * timeWarp;
                deltaV = fmin(fmax(deltaV, -0.5), 0.5);
                sat.velocity += v_direction * deltaV;
            }

            if (abs(targetInclination - currentInclination) > 0.5f)
            {
                Eigen::Vector3d r = sat.position;
                Eigen::Vector3d v = sat.velocity;
                Eigen::Vector3d normal = r.cross(v).normalized();
                double direction = (targetInclination - currentInclination > 0) ? 1.0f : -1.0f;
                double deltaV = direction * 0.08 * dt * timeWarp;
                sat.velocity += normal * deltaV;
            }

            if (t >= 1.0f)
            {
                isManeuvering = false;
            }
        }

        
        ImGui::Checkbox("Automatic Rotation", &auto_rotate);
        ImGui::Checkbox("Top-Down View", &use2p5D);

        ImGui::End();
        rlImGuiEnd();

        EndDrawing();
    }

    CloseWindow();
    return 0;
}

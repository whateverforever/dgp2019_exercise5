#include "viewer.h"
#include <iomanip>

using namespace surface_mesh;
typedef Surface_mesh Mesh;

bool isClose(float a, float b, float tol) {
    return std::abs(a - b) <= tol;
}

bool isClose(const Point &a, const Point &b, float tol=1e-2f) {
    return isClose(a[0], b[0], tol) && isClose(a[1], b[1], tol) && isClose(a[2], b[2], tol);
}

void print(const Point &p, bool brackets=true) {
    if (brackets) {
        std::cout << "[" << p[0] << ", " << p[1] << ", " << p[2] << "]" << std::endl;
    } else {
        std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
    }
}

int main(int argc, char **argv) {
    // cout.precision(15);

    //Loading the mesh.
    nanogui::init();
    nanogui::ref<Viewer> app = new Viewer(true);

    int n_successes = 0;
    int n_tests = 0;

    // Make sure we are in the bunny mesh
    app->loadMesh("../data/bunny.off");
    app->meshProcess();
    Surface_mesh &mesh = app->mesh;
    if (!(mesh.n_vertices() == 34835 && mesh.n_faces() == 69666)) {
        std::cout << "Wrong mesh was loaded." << std::endl;
        return 0;
    }


    // Get relevant buffers
    auto v_cste_weights_n = mesh.vertex_property<Point>("v:cste_weights_n");
    auto v_area_weights_n = mesh.vertex_property<Point>("v:area_weight_n");
    auto v_angle_weights_n = mesh.vertex_property<Point>("v:angle_weight_n");

    auto v_uniLaplace = mesh.vertex_property<Scalar>("v:uniLaplace", 0);
    auto v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    auto v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);

    auto vertex_valence = mesh.vertex_property<Scalar>("v:valence", 0);

    // Single out an interesting vertex
    Mesh::Vertex v_test;
    size_t k = 0;
    for (auto v: mesh.vertices()) {
        if (k == 101) {
            v_test = v;
        }
        k++;
    }

    Point n_test, n_ref, u_test, u_ref;

    std::cout << "======================================================================" << std::endl;
    std::cout << "Curvature test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "Part 1) Weighted vertex normals" << std::endl;

    std::cout << "  1.1: Constant weights ....................................... ";

    n_test = v_cste_weights_n[v_test];
    n_ref  = {-0.234669998288155f, -0.971998989582062f, 0.01216010004282f};
    if (isClose(n_test, n_ref)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(n_test);
    std::cout << "    reference: "; print(n_ref);
    n_tests++;

    std::cout << "  1.2: Area weights ........................................... ";

    n_test = v_area_weights_n[v_test];
    n_ref  = {0.0568871982395649f, -0.995039999485016f, -0.081599198281765f};
    if (isClose(n_test, n_ref)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(n_test);
    std::cout << "    reference: "; print(n_ref);
    n_tests++;

    std::cout << "  1.3: Angle weights .......................................... ";

    n_test = v_angle_weights_n[v_test];
    n_ref  = {-0.297583013772964f, -0.953437983989716f, -0.0489905998110771f};
    if (isClose(n_test, n_ref)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(n_test);
    std::cout << "    reference: "; print(n_ref);
    n_tests++;

    std::cout << std::endl;
    std::cout << "Part 2) Discrete curvature operators" << std::endl;

    std::cout << "  2.1: Uniform Laplace ........................................ ";

    u_test = {v_uniLaplace[v_test], app->min_uniLaplace, app->max_uniLaplace};
    u_ref  = {0.00417693983763456f, 4.10648993920404e-07f, 0.00598778016865253f};
    if (isClose(u_test, u_ref, 1e-3f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value/min/max:     "; print(u_test, false);
    std::cout << "    references:        "; print(u_ref, false);
    n_tests++;

    std::cout << "  2.2: Laplace-Beltrami ....................................... ";

    u_test = {v_curvature[v_test], app->min_mean_curvature, app->max_mean_curvature};
    u_ref  = {40.812801361084f, 0.0015746500575915f, 10226.7001953125f};

    if (isClose(u_test, u_ref, 1e-1f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value/min/max:     "; print(u_test, false);
    std::cout << "    references:        "; print(u_ref, false);
    n_tests++;

    std::cout << "  2.3: Gaussian ............................................... ";

    u_test = {v_gauss_curvature[v_test], app->min_gauss_curvature, app->max_gauss_curvature};
    u_ref  = {-5559.3701171875f, -1518350.f, 2178770.f};

    if (isClose(u_test, u_ref, 10.f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value/min/max:     "; print(u_test, false);
    std::cout << "    references:        "; print(u_ref, false);
    n_tests++;

    std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "PASSED " << n_successes << "/" << n_tests << " TESTS" << std::endl;

    nanogui::shutdown();
}

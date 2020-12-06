#include <timer.h>
#include <igl/get_seconds.h>

double timer(std::function<void(void)> func) {
    const double t_before = igl::get_seconds();
    func();
    const double t_after = igl::get_seconds();
    return t_after - t_before;
}
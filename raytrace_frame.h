#ifndef RAYTRACE_FRAME_H
#define RAYTRACE_FRAME_H

#include <vector>

bool raytrace_frame(const std::vector<double> &p,
                    const std::vector<double> &v,
                    double DX, double DY,
                    double &t_f,
                    std::vector<double> &C_f);

#endif
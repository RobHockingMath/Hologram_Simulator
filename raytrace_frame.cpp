#include "raytrace_frame.h"
#include <algorithm>
#include <cmath>
#include <cstddef>

// Helper functions for vector operations
static double norm(const std::vector<double> &a) {
    double s=0; for (auto val:a) s+=val*val;
    return std::sqrt(s);
}

static double dot(const std::vector<double> &a, const std::vector<double> &b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

static std::vector<double> normalize(const std::vector<double> &a) {
    double n=norm(a);
    if (n<1e-12) return {0,0,0};
    return {a[0]/n,a[1]/n,a[2]/n};
}

bool raytrace_frame(const std::vector<double> &p,
                    const std::vector<double> &v,
                    double DX, double DY,
                    double &t_f,
                    std::vector<double> &C_f)
{
    double dL = DX/40.0;

    std::vector<double> light_pos1 = {2*DX,2*DX,2*DY};
    std::vector<double> light_pos2 = {-2*DX,-2*DX,-2*DY};

    double h=2.0;

    std::vector<double> ts;
    std::vector<std::vector<double>> Ns;
    std::vector<std::vector<double>> ps;

    // Each intersection test similar to MATLAB code
    auto check_top = [&](double t1){
        if(t1>=0){
            double x1=p[1]+t1*v[1];
            double y1=p[2]+t1*v[2];
            if((( -dL<=x1 && x1<=0 && -dL<=y1 && y1<=DY+dL ) ||
                ( DX<=x1 && x1<=DX+dL && -dL<=y1 && y1<=DY+dL ) ||
                ( -dL<=y1 && y1<=0 && -dL<=x1 && x1<=DX+dL ) ||
                ( DY<=y1 && y1<=DY+dL && -dL<=x1 && x1<=DX+dL ))){
                std::vector<double> N={-1,0,0};
                std::vector<double> p_h={h,x1,y1};
                ts.push_back(t1);
                Ns.push_back(N);
                ps.push_back(p_h);
            }
        }
    };

    // t1 for top
    {
        double t1=(h-p[0])/v[0];
        check_top(t1);
    }

    auto check_side = [&](double t2, double zlow, double zhigh, double ylow, double yhigh,
                          const std::vector<double> &N_fixed, double px, double py) {
        // z1,x1,y1 meaning: in code z is p[0], x is p[1], y is p[2].
        // The MATLAB code variable naming is tricky, but we follow it directly.
        // Actually in code:
        // For vertical sides: we solve for x or y direction. We'll trust the original logic.
        if(t2>=0) {
            double z1 = p[0]+t2*v[0];
            double x1 = p[1]+t2*v[1];
            double y1 = p[2]+t2*v[2];
            // check if within bounds:
            if(z1>=zlow && z1<=zhigh && y1>=ylow && y1<=yhigh) {
                // intersection
                ts.push_back(t2);
                Ns.push_back(N_fixed);
                ps.push_back({z1,px,y1});
            }
        }
    };

    // Outer left side:
    {
        double t2=(-dL - p[1])/v[1];
        // (0<=z1 && z1<=h && -dL<=y1 && y1<=DY+dL)
        check_side(t2,0,h,-dL,DY+dL,{0,-1,0},-dL,0);
    }

    // Inner left side:
    {
        double t2=(0 - p[1])/v[1];
        // (0<=z1 && z1<=h && 0<=y1 && y1<=DY)
        check_side(t2,0,h,0,DY,{0,1,0},0,0);
    }

    // Outer right side:
    {
        double t2=(DX+dL - p[1])/v[1];
        // (0<=z1<=h && -dL<=y1<=DY+dL)
        check_side(t2,0,h,-dL,DY+dL,{0,1,0},DX+dL,0);
    }

    // Inner right side:
    {
        double t2=(DX - p[1])/v[1];
        // (0<=z1<=h && 0<=y1<=DY)
        check_side(t2,0,h,0,DY,{0,-1,0},DX,0);
    }

    // Outer bottom side:
    {
        double t2=(-dL - p[2])/v[2];
        // (0<=z1<=h && -dL<=x1<=DX+dL)
        // N=[0;0;-1], p_h=[z1;x1;-dL]
        if(t2>=0){
            double z1=p[0]+t2*v[0];
            double x1=p[1]+t2*v[1];
            double y1=p[2]+t2*v[2];
            if(z1>=0 && z1<=h && x1>=-dL && x1<=DX+dL) {
                ts.push_back(t2);
                Ns.push_back({0,0,-1});
                ps.push_back({z1,x1,-dL});
            }
        }
    }

    // Inner bottom side:
    {
        double t2=(0 - p[2])/v[2];
        // (0<=z1<=h && 0<=x1<=DX)
        if(t2>=0){
            double z1=p[0]+t2*v[0];
            double x1=p[1]+t2*v[1];
            double y1=p[2]+t2*v[2];
            if(z1>=0 && z1<=h && x1>=0 && x1<=DX) {
                ts.push_back(t2);
                Ns.push_back({0,0,1});
                ps.push_back({z1,x1,0});
            }
        }
    }

    // Outer top side:
    {
        double t2=(DY+dL - p[2])/v[2];
        // (0<=z1<=h && -dL<=x1<=DX+dL)
        if(t2>=0){
            double z1=p[0]+t2*v[0];
            double x1=p[1]+t2*v[1];
            double y1=p[2]+t2*v[2];
            if(z1>=0 && z1<=h && x1>=-dL && x1<=DX+dL){
                ts.push_back(t2);
                Ns.push_back({0,0,1});
                ps.push_back({z1,x1,DY+dL});
            }
        }
    }

    // Inner top side:
    {
        double t2=(DY - p[2])/v[2];
        // (0<=z1<=h && 0<=x1<=DX)
        if(t2>=0){
            double z1=p[0]+t2*v[0];
            double x1=p[1]+t2*v[1];
            double y1=p[2]+t2*v[2];
            if(z1>=0 && z1<=h && x1>=0 && x1<=DX) {
                ts.push_back(t2);
                Ns.push_back({0,0,-1});
                ps.push_back({z1,x1,DY});
            }
        }
    }

    if(!ts.empty()) {
        // find minimum t
        double t_min=ts[0]; size_t i_min=0;
        for (size_t i=1;i<ts.size();i++){
            if(ts[i]<t_min){ t_min=ts[i]; i_min=i;}
        }

        t_f=t_min;
        std::vector<double> N=Ns[i_min];
        std::vector<double> p_h=ps[i_min];

        // shading
        // delta1=light_pos1 - p_h
        std::vector<double> light_pos1={2*DX,2*DX,2*DY};
        std::vector<double> light_pos2={-2*DX,-2*DX,-2*DY};

        std::vector<double> delta1={light_pos1[0]-p_h[0],light_pos1[1]-p_h[1],light_pos1[2]-p_h[2]};
        std::vector<double> delta2={light_pos2[0]-p_h[0],light_pos2[1]-p_h[1],light_pos2[2]-p_h[2]};

        double shading1=std::max(dot(normalize(delta1),N),0.0);
        double shading2=std::max(dot(normalize(delta2),N),0.0);
        double shading=std::min(shading1+shading2,1.0);

        std::vector<double> C={255,255,255};
        C_f={shading*C[0],shading*C[1],shading*C[2]};
        return true;
    } else {
        t_f=-1;
        C_f={0,0,0};
        return false;
    }
}
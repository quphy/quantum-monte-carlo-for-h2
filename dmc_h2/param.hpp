#ifndef parameter
#define parameter

class Parameters{
public:
    inline static constexpr int num_ele = 2;
    inline static constexpr int num_nuc =2 ;
    inline static constexpr int niter =50000 ;
    inline static constexpr int eq_iter = 1000;
    inline static constexpr double dt = 0.01;
    inline static constexpr int max_walker =100000;
    inline static constexpr int avg_step = 1000 ;
    inline static constexpr int r_step = 1000 ;
    inline static constexpr double r_max =10.0;
    inline static constexpr int num_walker = 10000;
    inline static constexpr double R = 1.4;
};

#endif
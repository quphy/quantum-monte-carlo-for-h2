#ifndef parameter
#define parameter

class Parameters{
public:
    inline static constexpr int num_ele = 2;
    inline static constexpr int num_nuc =2 ;
    inline static constexpr double alpha = 2.0;
    inline static constexpr double gamma = 1.0;
    inline static constexpr int step = 100000;
    inline static constexpr int eq_step = 1000;
    inline static constexpr int num_walker = 1000;
    inline static constexpr int beta_step = 20;
    inline static constexpr double R = 1.4;
};

#endif
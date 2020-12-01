//
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
//

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <stdlib.h>
//#define PI 3.14159265

typedef boost::numeric::ublas::vector<double> vectord;
using namespace boost::numeric::odeint;
using namespace std;


struct k_ring
{
    double  m_omega;
    double  m_sigma;
    int     m_range;
    float m_alpha;
    int m_N;
    vectord m_coupling;

    k_ring( double omega = 2.0 , double sigma = 1.0 , int range = 1 , float alpha = M_PI-0.1 , int N = 5 ) :
    m_omega( omega ) , m_sigma ( sigma ) , m_range ( range ) ,
    m_alpha ( alpha ) , m_N ( N ) , m_coupling(N , 0.) {}

    void operator() (const vectord &x, vectord &dxdt, const double )
    {
        for( int i=0 ; i<m_N ; ++i ){
            //double a = i*2*M_PI/(m_N-1) - M_PI;
            m_coupling(i) = 0.;
            for( int j=0; j<m_N ; ++j ){
                //double b = j*2*M_PI/(m_N-1) - M_PI;
                float dist = abs(j-i);
                dist = abs(dist - round(dist/( (float) m_N ) ) * m_N);
                //m_coupling(i) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(j) - x(i) + m_alpha );
                //m_coupling(i) += (1/128)*exp(-(4/256)*dist)*sin( x(j) - x(i) + m_alpha );
                if(dist <= m_range && dist > 0){ //m_range
                    m_coupling(i) += (m_sigma/(2*m_range))*sin( x(j) - x(i) + m_alpha );
                }
            }
            dxdt(i) = m_omega + 2*M_PI*m_coupling(i)/(double) m_N;
        }
    }
};


struct observer
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    observer(ostream &out, const float &sim_time) :
        m_outfile(out) ,
        m_sim_time(sim_time){}
    
    void operator()(const vectord &x, double t) 
    {
        if(t > m_sim_time-10){
            for(size_t i=0; i < x.size(); ++i){
                if(fmod(t, 0.1) < 0.01){
                m_outfile << t << '\t' << i << '\t' << x(i) << '\n';
                }
            }
        }
    }
};
 
int main(int argc, char **argv)
{
    clock_t tStart = clock();
    
    //----------------------Reading parameter file-------------------------------//

        cout << "Reading parameter file..." << endl;
        string file = string(argv[1]); // "./parameters/" +

        boost::property_tree::ptree desc;
        boost::property_tree::json_parser::read_json(file, desc);

        const size_t   n     = desc.get<size_t>("n");
        const int      range_one   = desc.get<int>("range_one");
        const float    sig_one    = desc.get<float>("sig_one");
        const float    o_one    = desc.get<float>("o_one");
        const double   sim_time = desc.get<float>("sim_time");
        string         path    = desc.get<string>("path");
        
        cout << "Number of oscillators: "       << n   << '\n' <<
                "Coupling radius: "             << range_one    << '\n' <<
                "Coupling strength: "           << sig_one  << '\n' <<
                "Simulation time: "             << sim_time << '\n';

    //----------------------Setting initial conditions---------------------------//
    
    uniform_real_distribution<double> distribution(-0.5,0.5);
    random_device rd;
    default_random_engine generator(rd());
    
    //const size_t n = 256;
    //int range_one = 1;
    //float o_one = 1.3;
    //float sig_one = 1;
    float a = M_PI/2 - 0.18;
    //const double sim_time = 200.0;
    
    vectord state( n , 0. );
    
    for(size_t i = 0 ; i < n ; ++i )
    {
        double pos = i*2*M_PI/(n-1) - M_PI;
        double r = distribution(generator);
        state(i) = 6*r*exp(-0.76*pos*pos);
    }
    
    k_ring system(o_one, sig_one, range_one, a, n);

    string output_data = path;
    ofstream data_out(output_data);
    
    observer obs(data_out, sim_time);
    
    cout << "Starting integration..." << endl;
    
    adams_bashforth_moulton<5, vectord> abm;
    
    integrate_const(abm, system, state, 0.0, sim_time, 0.025, obs);

    data_out.close();

    cout << "Program successfully finished!" << '\n';
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

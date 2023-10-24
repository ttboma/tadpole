#include <topology.h>
#include <boost/graph/graphviz.hpp>

using CoordinateType        = double; 
using PointType			    = boost::polygon::point_data<CoordinateType>;
using VoronoiDiagramType    = boost::polygon::voronoi_diagram<CoordinateType>;

int main( int argc, char** argv ) 
{
    auto    ViaResults_FileName    = "via_info/500Via_Result.txt";
    auto    IoDrc_FileName         = "Io_drc/500Io_drc.txt";

    // parse and gnuplot input data 
    
    auto    pckg = syc::topology::v_03::package_data<>( ViaResults_FileName, IoDrc_FileName );
    auto    [ B, P, N ] =
            pckg.generate_fucked_up_case_using_bumps();
  
    for ( unsigned i = 0, i_end = B.size(); i != i_end; ++i )
    {
        boost::geometry::set<0>( B[i], 10 * boost::geometry::get<0>( B[i] ) ); 
        boost::geometry::set<1>( B[i], 10 * boost::geometry::get<1>( B[i] ) ); 
    }
    for ( unsigned i = 0, i_end = P.size(); i != i_end; ++i )
    {
        boost::geometry::set<0>( P[i], 10 * boost::geometry::get<0>( P[i] ) ); 
        boost::geometry::set<1>( P[i], 10 * boost::geometry::get<1>( P[i] ) ); 
    }
  
    auto    t1     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    sol    = syc::topology::v_03::planar_orderd_forest_solution<>( B, P, N );

    auto    t2     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

//  sol.find_MST_by_comparable_length();
//  sol.find_MST_by_A_star_heu_distance();
//  sol.find_MinLevelST();
    sol.find_MST_by_heuristic_method();
//  sol.find_MST_by_target_dist_from_root();

    auto    t3     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    cf     = syc::topology::v_03::circular_frame( sol );		
    
//  cf.the_basic_routing_algorithm();
//  cf.t_escape(); 
//  cf.t_escape2(); 
//  cf.t_escape3(); 
    cf.t_escape4(); 

    auto    t4     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    gp     = cf.realization( pckg.spacing * 30 );

    auto    t5     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    n      = cf.get_nets();

    auto    t6     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    double  T1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    double  T2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    double  T3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
    double  T4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();

    sol.gnuplot_svg( gp, "~/Desktop/xxx500", 1200, 1000 );

    std::cout << "# Time (microsecond)\n";
    std::cout << "Triangulation   = " << T1 << "\n";
    std::cout << "Find_MST        = " << T2 << "\n";
    std::cout << "Routing         = " << T3 << "\n";
    std::cout << "Realization     = " << T4 << "\n";
    std::cout << "\n"; 

    return 0;
}

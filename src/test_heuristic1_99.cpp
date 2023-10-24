#include <topology.h>

int main( int argc, char** argv ) 
{
    auto ViaResults_FileName    = "via_info/99Via_Result.txt";
    auto IoDrc_FileName         = "Io_drc/99Io_drc.txt";

    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( ViaResults_FileName, IoDrc_FileName );

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run the_basic_routing_algorithm and realize it.
     
    auto    t1      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    sol1    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    auto    t2      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    sol1.find_MST_by_comparable_length();
		
    auto    t3      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    cf1     = circular_frame( sol1 );		
	auto    cnt1    = cf1.the_basic_routing_algorithm();		

    auto    t4      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    gp1     = cf1.realization();

    auto    t5      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    n1      = cf1.get_nets();

    sol1.gnuplot( gp1 );

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run the_basic_routing_algorithm and realize it.
     
    auto    t6      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    sol2    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    auto    t7      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    sol2.find_MST_by_heuristic_method();
		
    auto    t8      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    cf2     = circular_frame( sol2 );		
	auto    cnt2    = cf2.the_basic_routing_algorithm();		

    auto    t9      = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    gp2     = cf2.realization();

    auto    t10     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    n2      = cf2.get_nets();

    sol2.gnuplot( gp2 );

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run t_escape2 and realize it.
     
    auto    t11     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    sol3    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    auto    t12     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    sol3.find_MST_by_comparable_length();

    auto    t13     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    cf3     = circular_frame( sol3 );		
	auto    cnt3    = cf3.t_escape2();		

    auto    t14     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    gp3     = cf3.realization();

    auto    t15     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    n3      = cf3.get_nets();

    sol3.gnuplot( gp3 );

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run t_escape2 and realize it.
     
    auto    t16     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    sol4    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    auto    t17     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    sol4.find_MST_by_heuristic_method();

    auto    t18     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    cf4     = circular_frame( sol4 );		
	auto    cnt4    = cf4.t_escape();		
    
    auto    t19     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    gp4     = cf4.realization();

    auto    t20     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );
    auto    n4      = cf4.get_nets();

    sol4.gnuplot( gp4 );
    sol4.gnuplot( gp4, "99.gif" );

    // output comparision
    
    auto    comp    = comparision( n3, gp3, n4, gp4 );
    comp.write_to_xlsx( "comparision_99" );

    std::cout << "triangulation:                 " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "(microsecond)\n";
    std::cout << "find_MST_by_comparable_length: " << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "(microsecond)\n";
    std::cout << "basic algo:                    " << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << "(microsecond)\n";
    std::cout << "realization:                   " << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() << "(microsecond)\n";
    std::cout << std::endl;

    std::cout << "triangulation:                 " << std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count() << "(microsecond)\n";
    std::cout << "find_MST_by_heuristic_method:  " << std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count() << "(microsecond)\n";
    std::cout << "basic algo:                    " << std::chrono::duration_cast<std::chrono::microseconds>(t9 - t8).count() << "(microsecond)\n";
    std::cout << "realization:                   " << std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count() << "(microsecond)\n";
    std::cout << std::endl;
    
    std::cout << "triangulation:                 " << std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count() << "(microsecond)\n";
    std::cout << "find_MST_by_comparable_length: " << std::chrono::duration_cast<std::chrono::microseconds>(t13 - t12).count() << "(microsecond)\n";
    std::cout << "t_escape:                      " << std::chrono::duration_cast<std::chrono::microseconds>(t14 - t13).count() << "(microsecond)\n";
    std::cout << "realization:                   " << std::chrono::duration_cast<std::chrono::microseconds>(t15 - t14).count() << "(microsecond)\n";
    std::cout << std::endl;

    std::cout << "triangulation:                 " << std::chrono::duration_cast<std::chrono::microseconds>(t17 - t16).count() << "(microsecond)\n";
    std::cout << "find_MST_by_heuristic_method:  " << std::chrono::duration_cast<std::chrono::microseconds>(t18 - t17).count() << "(microsecond)\n";
    std::cout << "t_escape:                      " << std::chrono::duration_cast<std::chrono::microseconds>(t19 - t18).count() << "(microsecond)\n";
    std::cout << "realization:                   " << std::chrono::duration_cast<std::chrono::microseconds>(t20 - t19).count() << "(microsecond)\n";
    std::cout << std::endl;

    return 0;
}

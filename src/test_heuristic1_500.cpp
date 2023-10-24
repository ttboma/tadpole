#include <topology.h>

template < typename NetSeq >
void sort_by_net_name_529 ( NetSeq& nets )
{
    std::sort( nets.begin(), nets.end(),
        [&] ( auto const& a, auto const& b )
        {
            auto    reg = std::regex( "(529_NET)([0-9]+)" );
            auto    sm_a  = std::smatch();
            auto    sm_b  = std::smatch();
            std::regex_search( a.name, sm_a, reg );
            std::regex_search( b.name, sm_b, reg );
            return std::stod( sm_a[2] ) < std::stod( sm_b[2] );
        }
    ); 
}

int main( int argc, char** argv ) 
{
    auto ViaResults_FileName    = "via_info/500Via_Result.txt";
    auto IoDrc_FileName         = "Io_drc/500Io_drc.txt";

    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( ViaResults_FileName, IoDrc_FileName );

    sort_by_net_name_529( pckg.nets );  // specifically needed for 500io 

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run the_basic_routing_algorithm and realize it.
     
    auto    sol2    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    sol2.find_MST_by_heuristic_method();
		
    auto    cf2     = circular_frame( sol2 );		
	auto    cnt2    = cf2.t_escape2();		
    auto    gp2     = cf2.realization();
    auto    n2      = cf2.get_nets();

    sol2.gnuplot( gp2 );

    return 0;
}

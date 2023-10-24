#include <topology.h>

auto generate_test_case( void )
{
    using namespace syc::topology::v_03;
	using point_type = typename planar_orderd_forest_solution< double >::point_type;

	std::vector< point_type > roots = 
	{
		point_type( 100, 100 ), point_type( 0, 100 ), point_type( -100, 100 ), point_type( -100, 0 ), 
		point_type( -100, -100 ), point_type( 100, -100 ) 	
    };
	std::vector< point_type > leafs = 
	{
        point_type( 20, 30 ), point_type( 30, 20 ), point_type( -10, 25 ), 
		point_type( -20, -10 ), point_type( -15, -40 )
	};
	std::vector< std::array< unsigned, 2 > > nets = 
	{
		{ 10, 8 }, { 9, 6 }
	};

	auto sol = planar_orderd_forest_solution< double >( roots, leafs, nets );

	return sol;
}
int main( int argc, char **argv )
{
    auto sol = generate_test_case();

    auto v_mom_mp = boost::get( vertex_mother, sol );
    v_mom_mp[6] = 1;
    v_mom_mp[7] = 6;
    v_mom_mp[8] = 6;
    v_mom_mp[9] = 3;
    v_mom_mp[10] = 3;

    auto e_idx_mp = boost::get( boost::edge_index, sol ); 
    auto const num_roots = 6;
    for ( auto [ vi, vi_end ] = vertices( sol ); vi != vi_end; ++vi )
    {
        if ( *vi < num_roots )          
            v_mom_mp[ *vi ] = *vi; 
        if ( v_mom_mp[ *vi ] != *vi )    
            e_idx_mp[ edge( *vi, v_mom_mp[ *vi ], sol ).first ] = 1; 
    }
    
    sol.gnuplot();

    auto cf = syc::topology::v_03::circular_frame( sol );		
    auto x1 = cf.make_slice(21,18);
    auto x2 = cf.make_slice(27,9);
    cf.txt_vitualizer();
    auto gp = cf.realization();
    sol.gnuplot(gp); // bug?? not working!!
//  sol.gnuplot_svg( gp, "~/Desktop/demo_6", 1200, 1000 );

    return 0;
}

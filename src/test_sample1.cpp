#include <topology.h>

// generate a small sample test case

auto generate_test_case( void )
{
    using namespace syc::topology::v_03;
	using point_type = typename planar_orderd_forest_solution< double >::point_type;

	std::vector< point_type > roots = 
	{
		point_type( 900, 900 ), point_type( 800, 900 ), point_type( 700, 900 ), point_type( 600, 900 ), 
		point_type( 500, 900 ), point_type( 400, 900 ), point_type( 300, 900 ), point_type( 200, 900 ), 
		point_type( 100, 900 ), point_type( 0, 900 ), point_type( 0, 0 ), point_type( 900, 0 ), 
		point_type( 900, 700 ), point_type( 900, 800 )
	};
	std::vector< point_type > leafs = 
	{
		point_type( 800, 600 ), point_type( 300, 600 ), point_type( 400, 450 ), point_type( 550, 350 ), 
		point_type( 700, 250 ), point_type( 250, 350 ), point_type( 250, 200 ), point_type( 400, 150 ),
		point_type( 50, 450 ), point_type( 850, 800 )
	};
	std::vector< std::array< unsigned, 2 > > nets = 
	{
		{ 17, 1 }, { 16, 7 }, { 15, 5 }, { 18, 4 }, { 19, 8 }, { 21, 6 }, { 20, 2 }, { 22, 13 }
	};

	auto sol = planar_orderd_forest_solution< double >( roots, leafs, nets );

	return sol;
}
int main( int argc, char** argv ) 
{
    auto sol1 = generate_test_case();

    // sol1.find_MST_by_comparable_length();
    sol1.find_MST_by_heuristic_method();
    
    sol1.print();
    sol1.gnuplot();
		
    auto    cf1     = syc::topology::v_02::circular_frame( sol1 );		
	auto    cnt1    = cf1.t_escape2();		
    auto    gp1     = cf1.realization();
    auto    n1      = cf1.get_nets();

    sol1.gnuplot( gp1 );
}

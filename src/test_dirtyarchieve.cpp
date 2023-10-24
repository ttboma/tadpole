#include <topology.h>
#include <simulated_anealing.h>

#include <random>
#include <ctime>

#include <boost/graph/graphviz.hpp>

void topology_test_1( void )
{
	// setting
	using coordinate_type			= double;
	using sketch_type				= syc::topology::model::v_01::topology_sketch_type<coordinate_type>;
	using point_type				= sketch_type::point_type;
	using sketchable_forest_type	= syc::topology::model::v_01::sketchable_forest<coordinate_type>;
	constexpr coordinate_type	packagesizeX		= 7000.000;
	constexpr coordinate_type	packagesizeY		= 7000.000;
	constexpr unsigned			num_nets_per_side	= 4;
	constexpr unsigned			num_nets			= num_nets_per_side * 4;
	// generating random points in a topology sketch TS
	constexpr coordinate_type spacingX = packagesizeX / ( num_nets_per_side );
	constexpr coordinate_type spacingY = packagesizeY / ( num_nets_per_side );
	std::uniform_real_distribution<coordinate_type> ux( -packagesizeX/2, packagesizeX/2 );
	std::uniform_real_distribution<coordinate_type> uy( -packagesizeY/2, packagesizeY/2 );
	std::default_random_engine e;
	std::vector<point_type> roots, leafs;
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2, packagesizeY/2 - i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2 + i * spacingX, -packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2, -packagesizeY/2 + i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2 - i * spacingX, packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets; ++i ) {
		leafs.emplace_back( ux( e ), uy( e ) );
	}
	// generating nets
	std::vector<unsigned> rand( num_nets );
		std::iota( rand.begin(), rand.end(), num_nets );
		std::shuffle( rand.begin(), rand.end(), std::mt19937{e()} );
	std::vector<std::pair<unsigned, unsigned>> nets;
	for ( unsigned i = 0; i != num_nets; ++i ) {
		nets.emplace_back( rand[ i ], i );
	}
	// topology routing 	
	// 1. construct a topology sketch, TS, with roots and leafs
	sketch_type TS( std::move( roots ), std::move( leafs ) );									
	// 2. perform t-escape algorithm with each of the nets whose source vertices are leafs and target vertices are roots 
	unsigned num_fail_nets = TS.t_escape( nets.begin(), nets.end() );		
	// 3, realize traces from the results
	auto traces = TS.realization();														 
	// 4. plot the results
	TS.plot( "results", traces );																
	// 5. analyze the results from realization
	for ( unsigned i = 0; i != traces.size(); ++i ) 
	{
		auto l = ( traces[i].empty() ) ? 0 : boost::geometry::length( traces[ i ] );
		std::cout << "length of net " << i << ": " << l << std::endl;
	}
}
void topology_test_2( void )
{
	// setting
	using coordinate_type			= double;
	using sketch_type				= syc::topology::model::v_01::topology_sketch_type<coordinate_type>;
	using point_type				= sketch_type::point_type;
	using sketchable_forest_type	= syc::topology::model::v_01::sketchable_forest<coordinate_type>;
	constexpr coordinate_type	packagesizeX		= 7000.000;
	constexpr coordinate_type	packagesizeY		= 7000.000;
	constexpr unsigned			num_nets_per_side	= 4;
	constexpr unsigned			num_nets			= num_nets_per_side * 4;
	// generating random points in a topology sketch TS
	constexpr coordinate_type spacingX = packagesizeX / ( num_nets_per_side );
	constexpr coordinate_type spacingY = packagesizeY / ( num_nets_per_side );
	std::uniform_real_distribution<coordinate_type> ux( -packagesizeX/2, packagesizeX/2 );
	std::uniform_real_distribution<coordinate_type> uy( -packagesizeY/2, packagesizeY/2 );
	std::default_random_engine e;
	std::vector<point_type> roots, leafs;
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2, packagesizeY/2 - i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2 + i * spacingX, -packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2, -packagesizeY/2 + i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2 - i * spacingX, packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets; ++i ) {
		leafs.emplace_back( ux( e ), uy( e ) );
	}

	// generating nets
	std::vector<unsigned> rand( num_nets );
		std::iota( rand.begin(), rand.end(), num_nets );
		std::shuffle( rand.begin(), rand.end(), std::mt19937{e()} );
	std::vector<std::pair<unsigned, unsigned>> nets;
	for ( unsigned i = 0; i != num_nets; ++i ) {
		nets.emplace_back( i, rand[ i ] );
	}

	// topology routing 	
	// 1. construct a topology sketch, TS, with roots and leafs
	sketch_type TS( std::move( roots ), std::move( leafs ) );									
	// 2. perform t-escape algorithm with each of the nets whose source vertices are leafs and target vertices are roots 
	unsigned num_fail_nets = TS.the_basic_routing_algorithm( nets.begin(), nets.end() );		
	// 3, realize traces from the results
	auto traces = TS.realization();														 
	// 4. plot the results
	TS.plot( "results", traces );																
	// 5. analyze the results from realization
	for ( unsigned i = 0; i != traces.size(); ++i ) 
	{
		auto l = ( traces[i].empty() ) ? 0 : boost::geometry::length( traces[ i ] );
		std::cout << "length of net " << i << ": " << l << std::endl;
	}
}
void topology_test_3( void )
{
	// setting
	using coordinate_type			= double;
	using sketch_type				= syc::topology::model::v_01::topology_sketch_type<coordinate_type>;
	using point_type				= sketch_type::point_type;
	using sketchable_forest_type	= syc::topology::model::v_01::sketchable_forest<coordinate_type>;
	constexpr coordinate_type	packagesizeX		= 7000.000;
	constexpr coordinate_type	packagesizeY		= 7000.000;
	constexpr unsigned			num_nets_per_side	= 10;
	std::default_random_engine	e;
	// generating random points in a topology sketch TS
	constexpr unsigned			num_nets			= num_nets_per_side * 4;
	constexpr coordinate_type spacingX = packagesizeX / ( num_nets_per_side );
	constexpr coordinate_type spacingY = packagesizeY / ( num_nets_per_side );
	std::uniform_real_distribution<coordinate_type> ux( -packagesizeX/2, packagesizeX/2 );
	std::uniform_real_distribution<coordinate_type> uy( -packagesizeY/2, packagesizeY/2 );
	std::vector<point_type> roots, leafs;
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2, packagesizeY/2 - i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2 + i * spacingX, -packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2, -packagesizeY/2 + i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2 - i * spacingX, packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets; ++i ) {
		leafs.emplace_back( ux( e ), uy( e ) );
	}

	// generating nets
//	std::vector<unsigned> rand( num_nets );
//		std::iota( rand.begin(), rand.end(), num_nets );
//		std::shuffle( rand.begin(), rand.end(), std::mt19937{e()} );
//	std::vector<std::pair<unsigned, unsigned>> nets;
//	for ( unsigned i = 0; i != num_nets; ++i ) {
//		nets.emplace_back( i, rand[ i ] );
//	}
	// topology routing 	
	// 1. construct a topology sketch, TS, with roots and leafs
	sketch_type TS( std::move( roots ), std::move( leafs ) );									
	// 2. perform t-escape algorithm with each of the nets whose source vertices are leafs and target vertices are roots 
	//unsigned num_fail_nets = TS.the_basic_routing_algorithm( nets.begin(), nets.end() );		
	// 3, realize traces from the results
	auto traces = TS.realization();														 
	// 4. plot the results
	TS.plot( "results", traces );																
	// 5. analyze the results from realization
	for ( unsigned i = 0; i != traces.size(); ++i ) 
	{
		auto l = ( traces[i].empty() ) ? 0 : boost::geometry::length( traces[ i ] );
		std::cout << "length of net " << i << ": " << l << std::endl;
	}
}

auto randomly_generate_a_case( double packagesizeX, double packagesizeY, unsigned num_nets_per_side ) 
{
	using Sol			= syc::topology::v_02::planar_orderd_forest_solution<double>;
	using point_type	= typename Sol::point_type;
	using return_type	= std::tuple<std::vector<point_type>, std::vector<point_type>, std::vector<std::array<std::size_t, 2>>>;

	return_type r;

	auto& roots = std::get<0>( r );
	auto& leafs = std::get<1>( r );

	double		spacingX = packagesizeX / ( num_nets_per_side );
	double		spacingY = packagesizeY / ( num_nets_per_side );
	unsigned	num_nets = num_nets_per_side * 4;

	std::uniform_real_distribution<double> ux( -packagesizeX/2, packagesizeX/2 );
	std::uniform_real_distribution<double> uy( -packagesizeY/2, packagesizeY/2 );

	std::default_random_engine e;

	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2, packagesizeY/2 - i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( -packagesizeX/2 + i * spacingX, -packagesizeY/2 );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2, -packagesizeY/2 + i * spacingX );
	}
	for ( unsigned i = 0; i != num_nets_per_side; ++i ) {
		roots.emplace_back( packagesizeX/2 - i * spacingX, packagesizeY/2 );
	}

	for ( unsigned i = 0; i != num_nets; ++i ) {
		leafs.emplace_back( ux( e ), uy( e ) );
	}
	
	auto& nets = std::get<2>( r );

	std::vector< std::size_t > temp( num_nets );

	std::iota( temp.begin(), temp.end(), num_nets );

	std::shuffle( temp.begin(), temp.end(), std::mt19937{e()} );

	for ( std::size_t i = 0; i != temp.size(); ++i ) 
	{
		nets.push_back( std::array<std::size_t, 2>{ i, temp[ i ] } );
	}

	return r;
}
void topology_test_4( void ) 
{
	auto [roots, leafs, nets] = 
		randomly_generate_a_case( 7000, 7000, 4 );	// produced 32 points 

	syc::topology::v_02::planar_orderd_forest_solution<double> 
		sol( roots, leafs, nets );

	auto sol2 = sol;

	std::default_random_engine 
		engine( 20 );
	std::uniform_int_distribution<std::size_t> 
		d( sol.num_roots, sol.num_leafs + sol.num_roots - 1 );
	for ( auto i = 0; i != 1000; ++i )
	{
		sol.switch_mother( d(engine) );
	}

	// sol
	sol.print();
	sol.gnuplot();
	std::cout << sol.cost() << std::endl;

	sol2.print();
	sol2.gnuplot();
	std::cout << sol2.cost() << std::endl;

}
void topology_test_5( void ) 
{
    auto [roots, leafs, nets] = 
      randomly_generate_a_case( 7000, 7000, 4 );	// produced 32 points 

    syc::topology::v_02::planar_orderd_forest_solution<double> 
      sol( roots, leafs, nets );

    auto thermal_equilibrium = [ i = 100 ]( auto currT ) mutable
    {
      if ( --i < 0 )
      {
        i = 100;
        return true;
      }
      else 
      {
        return false;
      }
    };
    auto decrease = [] ( auto& currT )
    {
      currT *= 0.85;
    };
    auto accept = 
      [ u		 = std::uniform_real_distribution< double >( 0, 1 ), 
        engine = std::default_random_engine(  std::chrono::system_clock::now().time_since_epoch().count() ) ] 
      ( auto currT, auto deltaC ) mutable 
    {
      return std::exp( - deltaC / currT ) > u( engine );
    };
    auto res = syc::SA::v_01::simulated_annealing(
      sol, 10000.0, 10.0, 
      thermal_equilibrium,
      decrease,
      accept
    );

      
    syc::topology::v_02::circular_frame	
          cf( res );		
    auto cnt    = cf.t_escape();		
    auto gp     = cf.realization();

    cf.txt_vitualizer();
    res.gnuplot( gp );
    std::cout << "cost: " << res.cost() << std::endl;

    auto print_results = [&] ( auto& cf, auto& gp, std::string words )
    {
        std::vector< std::tuple< unsigned, unsigned, double > > r;
        std::for_each( edges( cf.nets ).first, edges( cf.nets ).second, 
            [ &, i = int(0) ] ( auto e ) mutable
            {
                auto u = source( e, cf.nets );
                auto v = target( e, cf.nets );
                r.emplace_back( std::min( u , v ), std::max( u, v ), boost::geometry::length( gp[ i++ ] ) );
            }
        );
        std::sort( r.begin(), r.end(), 
            []( auto t1, auto t2 )
            {
                return std::min( std::get<0>( t1 ), std::get<1>( t1 ) ) < std::min( std::get<0>( t2 ), std::get<1>( t2 ) );
            } 
        );
        std::cout   << fmt::format( "{}\n", words )
                    << fmt::format( "{:-<{}}\n", "", words.size() );
        for ( auto const& result : r )
        {
            std::cout << fmt::format( "net [{:3},{:3}] of length: {:.2f}\n", std::get<0>(result), std::get<1>(result), std::get<2>(result) );
        }
        auto adv = std::accumulate( r.begin(), r.end(), 0, 
            [&] ( auto init, auto e )
            {
                return init + std::get<2>( e );
            }
        ) / r.size();
        std::cout << "adverage: " << adv << "\n";
        std::cout << std::endl;
    };
    print_results( cf, gp, "Results of the topology_escape algorithm" );
}
void topology_test_6( void ) 
{
	auto [roots, leafs, nets] = randomly_generate_a_case( 7000, 7000, 4 );	// produced 32 points 

	syc::topology::v_02::planar_orderd_forest_solution< double > 
        sol( roots, leafs, nets );
		
    // cost
    //
    std::cout << "\ncost: " << sol.cost() << "\n\n";

	syc::topology::v_02::circular_frame	
        cf( sol );		
	auto cnt = cf.t_escape2();		    //t_escape2 will show bugs: - -> + -> CW1
	auto gp = cf.realization();

    // output net and length 
    
    std::vector< std::tuple< unsigned, unsigned, double > > r;
    std::for_each( edges( cf.nets ).first, edges( cf.nets ).second, 
        [ &, i = int(0) ] ( auto e ) mutable
        {
            auto u = source( e, cf.nets );
            auto v = target( e, cf.nets );
            r.emplace_back( std::min( u , v ), std::max( u, v ), boost::geometry::length( gp[ i++ ] ) );
        }
    );
    std::sort( r.begin(), r.end(), 
        []( auto t1, auto t2 )
        {
            return std::min( std::get<0>( t1 ), std::get<1>( t1 ) ) < std::min( std::get<0>( t2 ), std::get<1>( t2 ) );
        } 
    );

	sol.gnuplot( gp );
    std::cout   << fmt::format( "Results of the basic routing algorithm\n" )
                << fmt::format( "--------------------------------------\n" );
    for ( auto const& result : r )
    {
        std::cout << fmt::format( "net [{:3},{:3}] of length: {:.2f}\n", std::get<0>(result), std::get<1>(result), std::get<2>(result) );
    }
    auto adv = std::accumulate( r.begin(), r.end(), 0, 
        [&] ( auto init, auto e )
        {
            return init + std::get<2>( e );
        }
    ) / r.size();
    std::cout << "adverage: " << adv << "\n";
   
    std::cout << std::endl;

}
void topology_test_6_2( void ) 
{
	auto [roots, leafs, nets]   = randomly_generate_a_case( 7000, 7000, 4 );	
	auto sol                    = syc::topology::v_02::planar_orderd_forest_solution< double >( roots, leafs, nets );
	auto cf                     = syc::topology::v_02::circular_frame( sol );		
//	auto cnt                    = cf.t_escape2();		                            // t_escape2 will show bugs: - -> + -> CW1
    auto    h = cf.make_slice( 2, 59 );
    cf.txt_vitualizer();
    cf.pretty_printer();
            h = cf.make_slice( h, 52 );
    cf.txt_vitualizer();
    cf.pretty_printer();
}
void topology_test_6_3( void ) 
{
    using namespace syc::topology::v_02;

    auto F = syc::topology::model::v_01::sketchable_forest< double >();
    add_edge( 4, 9, F );
    add_edge( 9, 6, F );
    add_edge( 6, 7, F );
    add_edge( 6, 8, F );
    add_edge( 8, 5, F );
    F.make_sketch();

    auto h1 = F.make_slice( 1, 19 );
    F.txt_vitualizer(); 
    F.pretty_printer(); 

    auto h2 = F.make_slice( h1, 12 );
    F.txt_vitualizer(); 
    F.pretty_printer(); 

    auto h3 = F.free_slice( h2 );
    F.txt_vitualizer(); 
    F.pretty_printer(); 

    auto h4 = F.free_slice( h3 );
    F.txt_vitualizer(); 
    F.pretty_printer(); 
}

auto generate_test_case1( void )
{
	using point_type = typename syc::topology::v_02::planar_orderd_forest_solution< double >::point_type;
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

	syc::topology::v_02::planar_orderd_forest_solution< double > sol( roots, leafs, nets );

	for ( auto i = 0; i != 1; ++i ) { sol.switch_mother( 15 ); }
	for ( auto i = 0; i != 2; ++i ) { sol.switch_mother( 16 ); }
	for ( auto i = 0; i != 5; ++i ) { sol.switch_mother( 22 ); }
	for ( auto i = 0; i != 2; ++i ) { sol.switch_mother( 19 ); }
	for ( auto i = 0; i != 1; ++i ) { sol.switch_mother( 21 ); }

	return sol;
}
void topology_test_7( void ) {
	// 這個case跑the_basic_routing_algorithm後，會在realization出問題
	// 且voronoi 轉triangulation有時會出現degenerate case，所以不會每個都是三角形。
	auto sol = generate_test_case1();

	syc::topology::v_02::circular_frame	cf1( sol );		
	auto cnt1 = cf1.t_escape();		
	auto gp1 = cf1.realization();
	sol.gnuplot( gp1 );
    std::vector< std::tuple< unsigned, unsigned, double > > r1;

	syc::topology::v_02::circular_frame	cf2( sol );		
	auto cnt2 = cf2.the_basic_routing_algorithm();		
	auto gp2 = cf2.realization();           
	sol.gnuplot( gp2 );
    std::vector< std::tuple< unsigned, unsigned, double > > r2;

    auto print_results = [&] ( auto& r, auto& cf, auto& gp, std::string words )
    {
        std::for_each( edges( cf.nets ).first, edges( cf.nets ).second, 
            [ &, i = int(0) ] ( auto e ) mutable
            {
                auto u = source( e, cf.nets );
                auto v = target( e, cf.nets );
                r.emplace_back( std::min( u , v ), std::max( u, v ), boost::geometry::length( gp[ i++ ] ) );
            }
        );
        std::sort( r.begin(), r.end(), 
            []( auto t1, auto t2 )
            {
                return std::min( std::get<0>( t1 ), std::get<1>( t1 ) ) < std::min( std::get<0>( t2 ), std::get<1>( t2 ) );
            } 
        );
        std::cout   << fmt::format( "{}\n", words )
                    << fmt::format( "{:-<{}}\n", "", words.size() );
        for ( auto const& result : r )
        {
            std::cout << fmt::format( "net [{:3},{:3}] of length: {:.2f}\n", std::get<0>(result), std::get<1>(result), std::get<2>(result) );
        }
        std::cout << std::endl;
    };
    print_results( r1, cf1, gp1, "Results of the topology escape routing algorithm" );
    print_results( r2, cf2, gp2, "Results of the basic routing algorithm" );

    std::cout << "sol cost: " << sol.cost() << std::endl;
}

syc::topology::v_02::planar_orderd_forest_solution<double> SA( syc::topology::v_02::planar_orderd_forest_solution<double>& sol )
{
    auto thermal_equilibrium = [ i = 100 ]( auto currT ) mutable
    {
      if ( --i < 0 )
      {
        i = 100;
        return true;
      }
      else 
      {
        return false;
      }
    };
    auto decrease = [] ( auto& currT )
    {
      currT *= 0.85;
    };
    auto accept = 
      [ u		 = std::uniform_real_distribution< double >( 0, 1 ), 
        engine = std::default_random_engine(  std::chrono::system_clock::now().time_since_epoch().count() ) ] 
      ( auto currT, auto deltaC ) mutable 
    {
      return std::exp( - deltaC / currT ) > u( engine );
    };
    return syc::SA::v_01::simulated_annealing( sol, 10000.0, 10.0, 
                                               thermal_equilibrium,
                                               decrease,
                                               accept                   );
}
void topology_test_8( std::string ViaResults_FileName, std::string IoDrc_FileName )     // ***
{
    using namespace syc::topology::v_02;

    // catch io keyword 
    auto reg    = std::regex( "(via_info/)([0-9]+)" );
    auto sm     = std::smatch();
                  std::regex_search( ViaResults_FileName, sm, reg );
    auto io     = std::string( sm[2] ); 

    // parse input data 
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    auto pckg = package_data<>( ViaResults_FileName, IoDrc_FileName );
    auto [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] = pckg.generate_case();
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
//  pckg.gnuplot();

    // run the basic routing algorithm and realize it
    auto sol    = planar_orderd_forest_solution<double>( roots, inner_bumps, inner_nets );
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    auto cf1    = circular_frame( sol );		
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
	auto cnt1   = cf1.the_basic_routing_algorithm();		
    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
    auto gp1    = cf1.realization();
    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
    sol.gnuplot( gp1 );
    std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();
    /*
    auto out1   = std::ofstream( io + "_basic.txt" );
    cf1.txt_vitualizer2( out1 );
    out1 << fmt::format( "cost: {}\n", sol.cost() );
    out1.close();
    */


    // run the topology escape routing algorithm and realize it
    auto sol2   = planar_orderd_forest_solution<double>( roots, inner_bumps, inner_nets );
    std::chrono::steady_clock::time_point t8 = std::chrono::steady_clock::now();
    auto res    = SA( sol ); 
    std::chrono::steady_clock::time_point t9 = std::chrono::steady_clock::now();
    
    // fucking ugly!!
    auto pfile  = std::ofstream( io + "_perturb_file.txt" );
    res.write_perturb( pfile );
    pfile.close();
    auto perturb_sample = std::ifstream( io + "_perturb_file.txt" ); 
    for ( unsigned v; perturb_sample >> v; ) { sol2.switch_mother( v ); }
    perturb_sample.close();
    auto clean_command = fmt::format( "rm {}\n", io + "_perturb_file.txt" );
    system( clean_command.c_str() );

    std::chrono::steady_clock::time_point t10 = std::chrono::steady_clock::now();
    auto cf2    = circular_frame( sol2 );		// Don't know why res is different from sol2
    std::chrono::steady_clock::time_point t11 = std::chrono::steady_clock::now();
    auto cnt2   = cf2.t_escape2();		
    std::chrono::steady_clock::time_point t12 = std::chrono::steady_clock::now();
    auto gp2    = cf2.realization();
    std::chrono::steady_clock::time_point t13 = std::chrono::steady_clock::now();
    sol2.gnuplot( gp2 );
    std::chrono::steady_clock::time_point t14 = std::chrono::steady_clock::now();
    /*
    auto out2   = std::ofstream( io + "_t_escape.txt" );
    cf2.txt_vitualizer2( out2 );
    out2 << fmt::format( "cost: {}\n", sol2.cost() );
    out2.close();
    */

    /*
    // There must be a bug why res is different from sol2. cf3 is different from cf2.
    auto cf3    = circular_frame( res );		
    auto cnt3   = cf3.t_escape2();		
    auto gp3    = cf3.realization();

    auto out3   = std::ofstream( io + "_t_escape_2.txt" );
    cf3.txt_vitualizer2( out3 );
    out3 << fmt::format( "cost: {}\n", res.cost() );
    out3.close();

    res.gnuplot( gp3 );
    */

    // compare results
    auto out3 = std::ofstream( io + "_results_comparasion" );
    auto print_comparasion = [&] ( std::ostream& os, auto const& cf1, auto const& gp1, auto const& cf2, auto const& gp2 )
    {
        auto results = [&] ( auto const& cf, auto const& gp )
        {
            std::vector< std::tuple< unsigned, unsigned, double > > r;
            std::for_each( edges( cf.nets ).first, edges( cf.nets ).second, 
                [ &, i = int(0) ] ( auto e ) mutable
                {
                    auto u = source( e, cf.nets );
                    auto v = target( e, cf.nets );
                    r.emplace_back( std::min( u , v ), std::max( u, v ), boost::geometry::length( gp[ i++ ] ) );
                }
            );
            std::sort( r.begin(), r.end(), 
                []( auto t1, auto t2 )
                {
                    return std::min( std::get<0>( t1 ), std::get<1>( t1 ) ) < std::min( std::get<0>( t2 ), std::get<1>( t2 ) );
                } 
            );
            return r;
        };
        auto adverage = [&] ( auto r )
        {
            return (long double)std::accumulate( r.begin(), r.end(), 0, 
                [&] ( long double init, auto e )
                {
                    return init + (long double)std::get<2>( e );
                }
            ) / r.size();
        };

        auto r1     = results( cf1, gp1 );
        auto avg1   = adverage( r1 );
        auto r2     = results( cf2, gp2 );
        auto avg2   = adverage( r2 );

        unsigned width[] = { 10, 40, 40, 20 };

        os << "The basic routing results\n\n";
        cf1.txt_vitualizer2( os );
        os << "\n";

        os << "The topology escape routing results\n\n";
        cf2.txt_vitualizer2( os );
        os << "\n";

        os  << fmt::format( 
            "{name: ^{w0}} | {basic: ^{w1}} | {topo: ^{w2}} | {normal: ^{w3}} \n"
            , fmt::arg( "w0", width[0] ), fmt::arg( "w1", width[1] ), fmt::arg( "w2", width[2] ), fmt::arg( "w3", width[3] )
            , fmt::arg( "name"      , io                        )                      
            , fmt::arg( "basic"     , "the basic algo"          )        
            , fmt::arg( "topo"      , "topology escape algo"    )  
            , fmt::arg( "normal"    , "comparasion"             )  
        );
        os  << fmt::format( 
            "{name:=^{w0}} | {basic:=^{w1}} | {topo:=^{w2}} | {normal:=^{w3}} \n"
            , fmt::arg( "w0", width[0] ), fmt::arg( "w1", width[1] ), fmt::arg( "w2", width[2] ), fmt::arg( "w3", width[3] )
            , fmt::arg( "name"      , "" )                      
            , fmt::arg( "basic"     , "" )        
            , fmt::arg( "topo"      , "" )  
            , fmt::arg( "normal"    , "" )
        );
        for ( unsigned i = 0; i != r1.size(); ++i )
        {
            auto n          = fmt::format( "[{}, {}]", std::get<0>( r1[i] ), std::get<1>( r1[i] ) );
            auto length1    = fmt::format ( "{:.2f}", std::get<2>( r1[i] ) );
            auto length2    = fmt::format ( "{:.2f}", std::get<2>( r2[i] ) );
            auto length     = fmt::format ( "{:.2f}", std::get<2>( r2[i] ) / std::get<2>( r1[i] ) );
            os  << fmt::format( 
                "{net: <{w0}} | {basic_length: ^{w1}} | {topo_length: ^{w2}} | {normal: ^{w3}} \n"
                , fmt::arg( "w0", width[0] ), fmt::arg( "w1", width[1] ), fmt::arg( "w2", width[2] ), fmt::arg( "w3", width[3] )
                , fmt::arg( "net"           , n         )                      
                , fmt::arg( "basic_length"  , length1   )        
                , fmt::arg( "topo_length"   , length2   )  
                , fmt::arg( "normal"        , length    )  
            );
        }
        os  << fmt::format( 
            "{name:=^{w0}} | {basic:=^{w1}} | {topo:=^{w2}} | {normal:=^{w3}} \n"
            , fmt::arg( "w0", width[0] ), fmt::arg( "w1", width[1] ), fmt::arg( "w2", width[2] ), fmt::arg( "w3", width[3] )
            , fmt::arg( "name"      , "" )                      
            , fmt::arg( "basic"     , "" )        
            , fmt::arg( "topo"      , "" )  
            , fmt::arg( "normal"    , "" )
        );
        auto AVG1 = fmt::format( "{:.2f}", avg1 );
        auto AVG2 = fmt::format( "{:.2f}", avg2 );
        auto AVG = fmt::format( "{:.2f}", avg2 / avg1 );
        os  << fmt::format( 
            "{name: ^{w0}} | {basic: ^{w1}} | {topo: ^{w2}} | {normal: ^{w3}} \n"
            , fmt::arg( "w0", width[0] ), fmt::arg( "w1", width[1] ), fmt::arg( "w2", width[2] ), fmt::arg( "w3", width[3] )
            , fmt::arg( "name"      , "average" )                      
            , fmt::arg( "basic"     , AVG1      )        
            , fmt::arg( "topo"      , AVG2      )  
            , fmt::arg( "normal"    , AVG       )  
        );
    };
    print_comparasion( out3, cf1, gp1, cf2, gp2 );


    std::cout   << fmt::format( "Parsing data:                  {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() );
    std::cout   << fmt::format( "Construct Initial Solution:    {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() );
    std::cout   << fmt::format( "Construct circular frame 1:    {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() );
    std::cout   << fmt::format( "the basic routing algorithm:   {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() );
    std::cout   << fmt::format( "geometry transformation 1:     {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() );
    std::cout   << fmt::format( "plot result 1:                 {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count() );
    std::cout   << "\n";
    std::cout   << fmt::format( "find solution by SA:           {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t9 - t8).count() );
    std::cout   << fmt::format( "Construct circular frame 2:    {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t11 - t10).count() );
    std::cout   << fmt::format( "topology escape routing:       {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count() );
    std::cout   << fmt::format( "geometry transformation 2:     {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t13 - t12).count() );
    std::cout   << fmt::format( "plot result 2:                 {}µs.\n", std::chrono::duration_cast<std::chrono::microseconds>(t14 - t13).count() );

}
void topology_test_9( void ) 
{
    using namespace syc::topology::v_02;

    auto pckg = package_data<>( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );
    auto [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] = pckg.generate_case();

    auto sol    = planar_orderd_forest_solution<double>( roots, inner_bumps, inner_nets );

    // perturb ...

    auto engine = std::default_random_engine( 243 );
    auto d      = std::uniform_int_distribution< std::size_t >( sol.num_roots, sol.num_leafs + sol.num_roots - 1 );
    for ( unsigned i = 0; i != 100; ++i )
    {
        while ( !sol.switch_mother( d( engine ) ) );
    }

    auto out_graph = std::ofstream( "sol" );
    boost::write_graphviz( out_graph, sol );
    out_graph.close();

    auto cf     = circular_frame( sol );		
	auto cnt    = cf.t_escape2();		
    auto out    = std::ofstream( "59_t_escape" );
    cf.txt_vitualizer2( out );


    auto gp     = cf.realization();     // bug are here
    sol.gnuplot( gp );
    std::cout << "cost: " << sol.cost() << std::endl;
}
void topology_test_10( void ) 
{
    using namespace syc::topology::v_02;

    auto pckg = package_data<>( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );
    auto [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] = pckg.generate_case();

    auto sol    = planar_orderd_forest_solution<double>( roots, inner_bumps, inner_nets );

    // perturb ...

    auto perturb_sample = std::ifstream( "59_perturb_file.txt" ); 
    for ( unsigned v; perturb_sample >> v; ) 
    { 
        sol.switch_mother( v ); 
    }

    auto cf     = circular_frame( sol );		
	auto cnt    = cf.t_escape2();		
    auto gp     = cf.realization();     

    auto out    = std::ofstream( "59_t_escape" );
    cf.txt_vitualizer2( out );

    sol.gnuplot();
    sol.gnuplot( gp );
}
void topology_test_11( std::string ViaResults_FileName, std::string IoDrc_FileName )    // *** 
{
    using namespace syc::topology::v_03;

    // catch io keyword 
    
    auto    reg = std::regex( "(via_info/)([0-9]+)" );
    auto    sm  = std::smatch();

    std::regex_search( ViaResults_FileName, sm, reg );

    auto    io  = std::string( sm[2] ); 

    // parse input data 
    
    auto    pckg = 
            package_data<>( ViaResults_FileName, IoDrc_FileName );

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run the_basic_routing_algorithm and realize it.
     
    auto    sol1    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    sol1.find_MST_by_comparable_length();
		
    auto    cf1     = syc::topology::v_02::circular_frame( sol1 );		
	auto    cnt1    = cf1.the_basic_routing_algorithm();		
    auto    gp1     = cf1.realization();
    auto    n1      = cf1.get_nets();

    sol1.gnuplot( gp1 );
    
    // Find a solution by marking MST on a delaunay triangulation undirected graph.
    // Construct circular frame from the solution.
    // Run t_escape2 and realize it.
     
    auto    sol2    = planar_orderd_forest_solution<>( roots, inner_bumps, inner_nets );

    sol2.find_MST_by_heuristic_method();

    auto    cf2     = syc::topology::v_02::circular_frame( sol2 );		
	auto    cnt2    = cf2.t_escape2();		
    auto    gp2     = cf2.realization();
    auto    n2      = cf2.get_nets();

    sol2.gnuplot( gp2 );

    // output comparision
    
    auto comp = comparision( n1, gp1, n2, gp2 );
    comp.write_to_xlsx( std::string( "comparision_" ) + io );
}

int main( int argc, char** argv ) 
{
     // topology_test_8( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );
     // topology_test_8( "via_info/99Via_Result.txt", "Io_drc/99Io_drc.txt" );
     topology_test_7();
    

    return 0;
}

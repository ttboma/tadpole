#include <topology.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

// avg_length

template < typename VecLinestring >
double avg_length( VecLinestring const& gp )
{
    auto L = std::vector< double >(); 
    for ( auto const line : gp ) 
    {
        L.push_back( boost::geometry::length(line) / gp.size() );
    }
    return std::accumulate( L.begin(), L.end(), 0.0 ); 
}

// MST_Policy

struct find_trees_by_random
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_trees_by_random();
    }
};
struct find_MaxST_by_comparable_length
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MaxST_by_comparable_length();
    }
};
struct find_MST_by_comparable_length
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MST_by_comparable_length();
    }
};
struct find_MST_by_A_star_heu_distance
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MST_by_A_star_heu_distance();
    }
};
struct find_MinLevelST
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MinLevelST();
    }
};
struct find_MST_by_heuristic_method
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MST_by_heuristic_method();
    }
};
struct find_MST_by_target_dist_from_root
{
    template < typename Sol >
    void find_MST( Sol& _s )
    {
        _s.find_MST_by_target_dist_from_root();
    }
};

// Routing_Policy

struct the_basic_routing_algorithm 
{
    template < typename CircularFrame >
    void apply_routing_policy( CircularFrame& _cf )
    {
	    _cf.the_basic_routing_algorithm();		
    }
};
struct topology_escape_routing_1
{
    template < typename CircularFrame >
    void apply_routing_policy( CircularFrame& _cf )
    {
	    _cf.t_escape();		
    }
};
struct topology_escape_routing_2
{
    template < typename CircularFrame >
    void apply_routing_policy( CircularFrame& _cf )
    {
	    _cf.t_escape2();		
    }
};
struct topology_escape_routing_3
{
    template < typename CircularFrame >
    void apply_routing_policy( CircularFrame& _cf )
    {
	    _cf.t_escape3();		
    }
};
struct topology_escape_routing_4
{
    template < typename CircularFrame >
    void apply_routing_policy( CircularFrame& _cf )
    {
	    _cf.t_escape4();		
    }
};

// work_flow 

template < 
    typename SimpleRing, typename InteriorPoint, typename Net, 
    typename MST_Policy, typename Routing_Policy 
>
auto work_flow( 
    SimpleRing const& B, InteriorPoint const& P, Net const& N, 
    MST_Policy F, Routing_Policy R, 
    std::string name
)
{
    using namespace syc::topology::v_03;

    auto    t1     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    sol    = planar_orderd_forest_solution<>( B, P, N );

    auto    t2     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    F.find_MST( sol );

    auto    t3     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    cf     = circular_frame( sol );		
    R.apply_routing_policy( cf ); 

    auto    t4     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    gp     = cf.realization();

    auto    t5     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    auto    n      = cf.get_nets();

    auto    t6     = std::chrono::steady_clock::time_point( std::chrono::steady_clock::now() );

    double  T1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    double  T2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    double  T3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
    double  T4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();

    sol.gnuplot_png( gp, name, 1200, 1000 );

    std::cout << "# Time (microsecond)\n";
    std::cout << "Triangulation   = " << T1 << "\n";
    std::cout << "Find_MST        = " << T2 << "\n";
    std::cout << "Routing         = " << T3 << "\n";
    std::cout << "Realization     = " << T4 << "\n";
    std::cout << "\n"; 

    return std::pair( avg_length( gp ), (T1 + T2 + T3 + T4) / 4 / 1000000 );
}

// Test

// 加 workfow 的流程
// 1. write a new function to find a forest
// 2. add to MST_Policy
// 3. add a new workflows below 
// 4. resize L and T, pics by 5, if any reslut of a routing policy is missing, used "NoResult" in pics

void Test_59( void )
{
    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // experiment on all conbination of MST policy and algorithm version
    
    auto workbook_name  = "result59.xlsx";
    auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
    auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
    auto pics           = std::vector<std::string>
    { 
        "result59_1"  , "result59_2"  , "result59_3"  , "result59_4"  , "result59_5"  ,
        "result59_6"  , "result59_7"  , "result59_8"  , "result59_9"  , "result59_10" ,
        "result59_11" , "result59_12" , "result59_13" , "result59_14" , "result59_15" ,
        "result59_16" , "result59_17" , "result59_18" , "result59_19" , "result59_20" ,
        "result59_21" , "result59_22" , "result59_23" , "result59_24" , "result59_25" ,
    };
    auto L = std::vector<double>( 25, 0 );
    auto T = std::vector<double>( 25, 0 );

    std::tie( L[0], T[0] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), the_basic_routing_algorithm(), pics[0] );
    std::tie( L[1], T[1] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_1(),   pics[1] );
    std::tie( L[2], T[2] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_2(),   pics[2] );
    std::tie( L[3], T[3] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_3(),   pics[3] );
    std::tie( L[4], T[4] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_4(),   pics[4] );
                       
    std::tie( L[5], T[5] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), the_basic_routing_algorithm(),  pics[5] ); 
    std::tie( L[6], T[6] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_1(), pics[6] );
    std::tie( L[7], T[7] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_2(), pics[7] );
    std::tie( L[8], T[8] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_3(), pics[8] );
    std::tie( L[9], T[9] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_4(), pics[9] );

    std::tie( L[10], T[10] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), the_basic_routing_algorithm(),  pics[10] ); 
    std::tie( L[11], T[11] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_1(), pics[11] );
    std::tie( L[12], T[12] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_2(), pics[12] );
    std::tie( L[13], T[13] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_3(), pics[13] );
    std::tie( L[14], T[14] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_4(), pics[14] );

    std::tie( L[15], T[15] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), the_basic_routing_algorithm(),  pics[15] );
    std::tie( L[16], T[16] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_1(),    pics[16] );
    std::tie( L[17], T[17] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_2(),    pics[17] );
    std::tie( L[18], T[18] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_3(),    pics[18] );
    std::tie( L[19], T[19] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_4(),    pics[19] );

    std::tie( L[20], T[20] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), the_basic_routing_algorithm(),  pics[20] );
    std::tie( L[21], T[21] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_1(),    pics[21] );
    std::tie( L[22], T[22] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_2(),    pics[22] );
    std::tie( L[23], T[23] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_3(),    pics[23] );
    std::tie( L[24], T[24] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_4(),    pics[24] );

    auto out = std::ofstream( "./src/write_to_xlsx_argv.py" );
    fmt::print( out, "workbook_name = '{}'  \n", workbook_name );
    fmt::print( out, "Routing_Policy = {}   \n", Routing_Policy );
    fmt::print( out, "MST_Policy = {}       \n", MST_Policy );
    fmt::print( out, "L = {}                \n", L );
    fmt::print( out, "T = {}                \n", T );
    fmt::print( out, "P = {}                \n", pics );
    out.close();
    system( "python3 ./src/write_to_xlsx.py" );
    system( "rm ./src/write_to_xlsx_argv.py" );
}
void Test_99( void )
{
    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( "via_info/99Via_Result.txt", "Io_drc/99Io_drc.txt" );

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // experiment on all conbination of MST policy and algorithm version
    
    auto workbook_name  = "result99.xlsx";
    auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
    auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
    auto L = std::vector<double>( 25, 0 );
    auto T = std::vector<double>( 25, 0 );
    auto pics = std::vector<std::string>
    { 
        "result99_1"    , "result99_2"  , "result99_3"  , "result99_4"  , "result99_5"  ,
        "./src/NoResult", "result99_7"  , "result99_8"  , "result99_9"  , "result99_10" ,
        "result99_11"   , "result99_12" , "result99_13" , "result99_14" , "result99_15" ,
        "result99_16"   , "result99_17" , "result99_18" , "result99_19" , "result99_20" ,
        "result99_21"   , "result99_22" , "result99_23" , "result99_24" , "result99_25" ,
    };

    std::tie( L[0], T[0] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), the_basic_routing_algorithm(), pics[0] );
    std::tie( L[1], T[1] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_1(),   pics[1] );
    std::tie( L[2], T[2] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_2(),   pics[2] );
    std::tie( L[3], T[3] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_3(),   pics[3] );
    std::tie( L[4], T[4] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_4(),   pics[4] );
                       
//    std::tie( L[5], T[5] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), the_basic_routing_algorithm(),  pics[5] ); // seg fault
    std::tie( L[6], T[6] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_1(), pics[6] );
    std::tie( L[7], T[7] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_2(), pics[7] );
    std::tie( L[8], T[8] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_3(), pics[8] );
    std::tie( L[9], T[9] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_4(), pics[9] );

    std::tie( L[10], T[10] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), the_basic_routing_algorithm(),  pics[10] ); 
    std::tie( L[11], T[11] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_1(), pics[11] );
    std::tie( L[12], T[12] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_2(), pics[12] );
    std::tie( L[13], T[13] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_3(), pics[13] );
    std::tie( L[14], T[14] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_4(), pics[14] );

    std::tie( L[15], T[15] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), the_basic_routing_algorithm(),  pics[15] );
    std::tie( L[16], T[16] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_1(),    pics[16] );
    std::tie( L[17], T[17] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_2(),    pics[17] );
    std::tie( L[18], T[18] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_3(),    pics[18] );
    std::tie( L[19], T[19] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_4(),    pics[19] );

    std::tie( L[20], T[20] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), the_basic_routing_algorithm(),  pics[20] );
    std::tie( L[21], T[21] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_1(),    pics[21] );
    std::tie( L[22], T[22] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_2(),    pics[22] );
    std::tie( L[23], T[23] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_3(),    pics[23] );
    std::tie( L[24], T[24] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_4(),    pics[24] );

    auto out = std::ofstream( "./src/write_to_xlsx_argv.py" );
    fmt::print( out, "workbook_name = '{}'  \n", workbook_name );
    fmt::print( out, "Routing_Policy = {}   \n", Routing_Policy );
    fmt::print( out, "MST_Policy = {}       \n", MST_Policy );
    fmt::print( out, "L = {}                \n", L );
    fmt::print( out, "T = {}                \n", T );
    fmt::print( out, "P = {}                \n", pics );
    out.close();
    system( "python3 ./src/write_to_xlsx.py" );
    system( "rm ./src/write_to_xlsx_argv.py" );
}
void Test_256( void )
{
    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( "via_info/256Via_Result.txt", "Io_drc/256Io_drc.txt" );

    std::sort( pckg.nets.begin(), pckg.nets.end(),
        [&] ( auto const& a, auto const& b )
        {
            auto    reg = std::regex( "(NET_)(([0-9]*[.])?[0-9]+)" );
            auto    sm_a  = std::smatch();
            auto    sm_b  = std::smatch();
            std::regex_search( a.name, sm_a, reg );
            std::regex_search( b.name, sm_b, reg );
            return std::stod( sm_a[2] ) < std::stod( sm_b[2] );
        }
    ); 

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // experiment on all conbination of MST policy and algorithm version
    
    auto workbook_name  = "result256.xlsx";
    auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
    auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
    auto L = std::vector<double>( 25, 0 );
    auto T = std::vector<double>( 25, 0 );
    auto pics = std::vector<std::string>
    { 
        "./src/NoResult", "result256_2"  , "result256_3"  , "result256_4"  , "result256_5"  ,
        "result256_6"   , "result256_7"  , "result256_8"  , "result256_9"  , "result256_10" ,
        "./src/NoResult", "result256_12" , "result256_13" , "result256_14" , "result256_15" ,
        "result256_16"  , "result256_17" , "result256_18" , "result256_19" , "result256_20" ,
        "result256_21"  , "result256_22" , "result256_23" , "result256_24" , "result256_25" ,
    };

//    std::tie( L[0], T[0] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), the_basic_routing_algorithm(), pics[0] );
    std::tie( L[1], T[1] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_1(),   pics[1] );
    std::tie( L[2], T[2] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_2(),   pics[2] );
    std::tie( L[3], T[3] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_3(),   pics[3] );
    std::tie( L[4], T[4] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_comparable_length(), topology_escape_routing_4(),   pics[4] );
                       
    std::tie( L[5], T[5] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), the_basic_routing_algorithm(),  pics[5] ); 
    std::tie( L[6], T[6] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_1(), pics[6] );
    std::tie( L[7], T[7] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_2(), pics[7] );
    std::tie( L[8], T[8] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_3(), pics[8] );
    std::tie( L[9], T[9] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_A_star_heu_distance(), topology_escape_routing_4(), pics[9] );

//    std::tie( L[10], T[10] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), the_basic_routing_algorithm(),  pics[10] ); 
    std::tie( L[11], T[11] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_1(), pics[11] );
    std::tie( L[12], T[12] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_2(), pics[12] );
    std::tie( L[13], T[13] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_3(), pics[13] );
    std::tie( L[14], T[14] ) = work_flow( roots, inner_bumps, inner_nets, find_MinLevelST(), topology_escape_routing_4(), pics[14] );

    std::tie( L[15], T[15] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), the_basic_routing_algorithm(),  pics[15] );
    std::tie( L[16], T[16] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_1(),    pics[16] );
    std::tie( L[17], T[17] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_2(),    pics[17] );
    std::tie( L[18], T[18] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_3(),    pics[18] );
    std::tie( L[19], T[19] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_4(),    pics[19] );

    std::tie( L[20], T[20] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), the_basic_routing_algorithm(),  pics[20] );
    std::tie( L[21], T[21] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_1(),    pics[21] );
    std::tie( L[22], T[22] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_2(),    pics[22] );
    std::tie( L[23], T[23] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_3(),    pics[23] );
    std::tie( L[24], T[24] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_target_dist_from_root(), topology_escape_routing_4(),    pics[24] );

    auto out = std::ofstream( "./src/write_to_xlsx_argv.py" );
    fmt::print( out, "workbook_name = '{}'  \n", workbook_name );
    fmt::print( out, "Routing_Policy = {}   \n", Routing_Policy );
    fmt::print( out, "MST_Policy = {}       \n", MST_Policy );
    fmt::print( out, "L = {}                \n", L );
    fmt::print( out, "T = {}                \n", T );
    fmt::print( out, "P = {}                \n", pics );
    out.close();
    system( "python3 ./src/write_to_xlsx.py" );
    system( "rm ./src/write_to_xlsx_argv.py" );
}

void Random_Test_59( void )
{
    using namespace syc::topology::v_03;

    // parse input data 
    
    auto    pckg = 
            package_data<>( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );

    auto    [ roots, inner_bumps, inner_nets, outer_bumps, outer_nets ] =
            pckg.generate_case_using_bumps();

    // experiment on all conbination of MST policy and algorithm version
    
    auto workbook_name  = "result59.xlsx";
    auto Routing_Policy = std::vector<std::string>{ "Basic", "T_Escape1", "T_Escape2" , "T_Escape3" , "T_Escape4" };
    auto MST_Policy     = std::vector<std::string>{ "Rand1", "Rand2"    , "Rand3"    , "Rand4" , "Rand5"   };
    auto pics           = std::vector<std::string>
    { 
        "result59_1"  , "result59_2"  , "result59_3"  , "result59_4"  , "result59_5"  ,
        "result59_6"  , "result59_7"  , "result59_8"  , "result59_9"  , "result59_10" ,
        "result59_11" , "result59_12" , "result59_13" , "result59_14" , "result59_15" ,
        "result59_16" , "result59_17" , "result59_18" , "result59_19" , "result59_20" ,
        "result59_21" , "result59_22" , "result59_23" , "result59_24" , "result59_25" ,
    };
    auto L = std::vector<double>( 25, 0 );
    auto T = std::vector<double>( 25, 0 );

    std::tie( L[0], T[0] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), the_basic_routing_algorithm(), pics[0] );
    std::tie( L[1], T[1] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_1(),   pics[1] );
    std::tie( L[2], T[2] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_2(),   pics[2] );
    std::tie( L[3], T[3] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_3(),   pics[3] );
    std::tie( L[4], T[4] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_4(),   pics[4] );
                       
    std::tie( L[5], T[5] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), the_basic_routing_algorithm(),  pics[5] ); 
    std::tie( L[6], T[6] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_1(), pics[6] );
    std::tie( L[7], T[7] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_2(), pics[7] );
    std::tie( L[8], T[8] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_3(), pics[8] );
    std::tie( L[9], T[9] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_4(), pics[9] );

    std::tie( L[10], T[10] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), the_basic_routing_algorithm(),  pics[10] ); 
    std::tie( L[11], T[11] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_1(), pics[11] );
    std::tie( L[12], T[12] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_2(), pics[12] );
    std::tie( L[13], T[13] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_3(), pics[13] );
    std::tie( L[14], T[14] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_4(), pics[14] );

    std::tie( L[15], T[15] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), the_basic_routing_algorithm(),  pics[15] );
    std::tie( L[16], T[16] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_1(),    pics[16] );
    std::tie( L[17], T[17] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_2(),    pics[17] );
    std::tie( L[18], T[18] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_3(),    pics[18] );
    std::tie( L[19], T[19] ) = work_flow( roots, inner_bumps, inner_nets, find_trees_by_random(), topology_escape_routing_4(),    pics[19] );

    std::tie( L[20], T[20] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), the_basic_routing_algorithm(),  pics[20] );
    std::tie( L[21], T[21] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_1(),    pics[21] );
    std::tie( L[22], T[22] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_2(),    pics[22] );
    std::tie( L[23], T[23] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_3(),    pics[23] );
    std::tie( L[24], T[24] ) = work_flow( roots, inner_bumps, inner_nets, find_MST_by_heuristic_method(), topology_escape_routing_4(),    pics[24] );

    auto out = std::ofstream( "write_to_xlsx_argv.py" );
    fmt::print( out, "workbook_name = '{}'  \n", workbook_name );
    fmt::print( out, "Routing_Policy = {}   \n", Routing_Policy );
    fmt::print( out, "MST_Policy = {}       \n", MST_Policy );
    fmt::print( out, "L = {}                \n", L );
    fmt::print( out, "T = {}                \n", T );
    fmt::print( out, "P = {}                \n", pics );
    out.close();
    system( "python3 write_to_xlsx.py" );
    system( "rm write_to_xlsx_argv.py" );
}

int main( int argc, char** argv ) 
{
    Test_59();
    Test_99();
    Test_256();

    return 0;
}

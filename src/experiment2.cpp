#include <topology.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>


namespace MST_Policy
{
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
}

namespace Routing_Policy
{
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
}

namespace Test
{
    //  加 workfow 的流程:
    //  1. write a new function to find a forest
    //  2. add to MST_Policy
    //  3. add a new workflows below 
    //  4. resize L and T, pics by 5, if any reslut of a routing policy is missing, used "NoResult" in pics

    namespace detail
    {
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
        template < typename SimpleRing, typename InteriorPoint, typename Net, 
                   typename Spacing, typename MST_Policy, typename Routing_Policy >
        auto work_flow( SimpleRing const& B, InteriorPoint const& P, Net const& N,
                        Spacing s, MST_Policy F, Routing_Policy R, std::string name )
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

            auto    gp     = cf.realization( s );

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

            return std::pair( avg_length( gp ), (T1 + T2 + T3 + T4) / 4 );
        }
    }    
      
    void io59( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/59Via_Result.txt", "Io_drc/59Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

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

        std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
        std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
        std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
        std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
        std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
        std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

        std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
        std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
        std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

        std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

        std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
    void io99( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/99Via_Result.txt", "Io_drc/99Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

        // experiment on all conbination of MST policy and algorithm version
        
        auto workbook_name  = "result99.xlsx";
        auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
        auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
        auto pics           = std::vector<std::string>
        { 
            "./src/NoResult" , "result99_2"  , "result99_3"  , "result99_4"  , "result99_5"  ,
            "./src/NoResult" , "result99_7"  , "result99_8"  , "result99_9"  , "result99_10" ,
            "result99_11"    , "result99_12" , "result99_13" , "result99_14" , "result99_15" ,
            "result99_16"    , "result99_17" , "result99_18" , "result99_19" , "result99_20" ,
            "result99_21"    , "result99_22" , "result99_23" , "result99_24" , "result99_25" ,
        };
        auto L = std::vector<double>( 25, 0 );
        auto T = std::vector<double>( 25, 0 );

//      std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
        std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
        std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
//      std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
        std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
        std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

        std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
        std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
        std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

        std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

        std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
    void io256( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/256Via_Result.txt", "Io_drc/256Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

        // experiment on all conbination of MST policy and algorithm version
        
        auto workbook_name  = "result256.xlsx";
        auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
        auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
        auto pics           = std::vector<std::string>
        { 
            "./src/NoResult" , "result256_2"  , "result256_3"  , "result256_4"  , "result256_5"  ,
            "./src/NoResult" , "result256_7"  , "result256_8"  , "result256_9"  , "result256_10" ,
            "./src/NoResult" , "result256_12" , "result256_13" , "result256_14" , "result256_15" ,
            "result256_16"   , "result256_17" , "result256_18" , "result256_19" , "result256_20" ,
            "result256_21"   , "result256_22" , "result256_23" , "result256_24" , "result256_25" ,
        };
        auto L = std::vector<double>( 25, 0 );
        auto T = std::vector<double>( 25, 0 );

//      std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
        std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
        std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
//      std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
        std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
        std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

//      std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
        std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
        std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

        std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

        std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
    void io300( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/300Via_Result.txt", "Io_drc/300Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

        // experiment on all conbination of MST policy and algorithm version
        
        auto workbook_name  = "result300.xlsx";
        auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
        auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
        auto pics           = std::vector<std::string>
        { 
            "result300_1"  , "result300_2"  , "result300_3"  , "result300_4"  , "result300_5"  ,
            "result300_6"  , "result300_7"  , "result300_8"  , "result300_9"  , "result300_10" ,
            "result300_11" , "result300_12" , "result300_13" , "result300_14" , "result300_15" ,
            "result300_16" , "result300_17" , "result300_18" , "result300_19" , "result300_20" ,
            "result300_21" , "result300_22" , "result300_23" , "result300_24" , "result300_25" ,
        };
        auto L = std::vector<double>( 25, 0 );
        auto T = std::vector<double>( 25, 0 );

        std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
        std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
        std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
        std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
        std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
        std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

        std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
        std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
        std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

        std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

        std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
    void io400( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/400Via_Result.txt", "Io_drc/400Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

        // experiment on all conbination of MST policy and algorithm version
        // result400_6 takes too long to plot
        auto workbook_name  = "result400.xlsx";
        auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
        auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
        auto pics           = std::vector<std::string>
        { 
            "./src/NoResult" , "result400_2"  , "result400_3"  , "result400_4"  , "result400_5"  ,
            "./src/NoResult" , "result400_7"  , "result400_8"  , "result400_9"  , "result400_10" ,
            "./src/NoResult" , "result400_12" , "result400_13" , "result400_14" , "result400_15" ,
            "result400_16"   , "result400_17" , "result400_18" , "result400_19" , "result400_20" ,
            "result400_21"   , "result400_22" , "result400_23" , "result400_24" , "result400_25" ,
        };
        auto L = std::vector<double>( 25, 0 );
        auto T = std::vector<double>( 25, 0 );

//      std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
        std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
        std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
//      std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
        std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
        std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

//      std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
        std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
        std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

        std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

        std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
    void io500( void )
    {
        using namespace syc::topology::v_03;

        // parse input data 
        
        auto    pckg = 
                package_data<>( "via_info/500Via_Result.txt", "Io_drc/500Io_drc.txt" );

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
        auto    s = pckg.spacing * 30;

        // experiment on all conbination of MST policy and algorithm version
        // result500_6 might have a result
        
        auto workbook_name  = "result500.xlsx";
        auto Routing_Policy = std::vector<std::string>{ "Basic"   , "T_Escape1"     , "T_Escape2" , "T_Escape3" , "T_Escape4"        };
        auto MST_Policy     = std::vector<std::string>{ "Eud_dis" , "Astar_Eud_dis" , "MinLevel"  , "delta_t"   , "delta_t_from_r"   };
        auto pics           = std::vector<std::string>
        { 
            "./src/NoResult" , "./src/NoResult"  , "./src/NoResult"  , "result500_4"  , "result500_5"  ,
            "./src/NoResult" , "./src/NoResult"  , "./src/NoResult"  , "result500_9"  , "result500_10" ,
            "./src/NoResult" , "./src/NoResult" , "./src/NoResult" , "result500_14" , "result500_15" ,
            "./src/NoResult" , "result500_17" , "result500_18" , "result500_19" , "result500_20" ,
            "./src/NoResult" , "result500_22" , "result500_23" , "result500_24" , "result500_25" ,
        };
        auto L = std::vector<double>( 25, 0 );
        auto T = std::vector<double>( 25, 0 );

//      std::tie( L[0], T[0] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::the_basic_routing_algorithm(), pics[0] );
//      std::tie( L[1], T[1] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_1(),   pics[1] );
//      std::tie( L[2], T[2] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_2(),   pics[2] );
        std::tie( L[3], T[3] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_3(),   pics[3] );
        std::tie( L[4], T[4] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_comparable_length(), Routing_Policy::topology_escape_routing_4(),   pics[4] );
                           
//      std::tie( L[5], T[5] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::the_basic_routing_algorithm(),  pics[5] ); 
//      std::tie( L[6], T[6] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_1(), pics[6] );
//      std::tie( L[7], T[7] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_2(), pics[7] );
        std::tie( L[8], T[8] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_3(), pics[8] );
        std::tie( L[9], T[9] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_A_star_heu_distance(), Routing_Policy::topology_escape_routing_4(), pics[9] );

//      std::tie( L[10], T[10] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::the_basic_routing_algorithm(),  pics[10] ); 
//      std::tie( L[11], T[11] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_1(), pics[11] );
//      std::tie( L[12], T[12] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_2(), pics[12] );
        std::tie( L[13], T[13] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_3(), pics[13] );
        std::tie( L[14], T[14] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MinLevelST(), Routing_Policy::topology_escape_routing_4(), pics[14] );

//      std::tie( L[15], T[15] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::the_basic_routing_algorithm(),  pics[15] );
        std::tie( L[16], T[16] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_1(),    pics[16] );
        std::tie( L[17], T[17] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_2(),    pics[17] );
        std::tie( L[18], T[18] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_3(),    pics[18] );
        std::tie( L[19], T[19] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_heuristic_method(), Routing_Policy::topology_escape_routing_4(),    pics[19] );

//      std::tie( L[20], T[20] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::the_basic_routing_algorithm(),  pics[20] );
        std::tie( L[21], T[21] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_1(),    pics[21] );
        std::tie( L[22], T[22] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_2(),    pics[22] );
        std::tie( L[23], T[23] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_3(),    pics[23] );
        std::tie( L[24], T[24] ) = detail::work_flow( B, P, N, s, MST_Policy::find_MST_by_target_dist_from_root(), Routing_Policy::topology_escape_routing_4(),    pics[24] );

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
}

int main( int argc, char** argv ) 
{
    Test::io59();
    Test::io99();
    Test::io256();
    Test::io300();
    Test::io400();
    Test::io500();
    
    return 0;
}

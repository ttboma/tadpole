#ifndef _SYC_TOPOLOGY_H
#define _SYC_TOPOLOGY_H

// standard library deps
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <tuple>
#include <map>
#include <memory>
#include <utility>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <regex>
// boost graph library deps
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/property.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp> 
enum vertex_diameter_t		{ vertex_diameter 		};
enum vertex_width_t			{ vertex_width 			};
enum vertex_coordinate_t	{ vertex_coordinate 	};
enum vertex_topoInstance_t	{ vertex_topoInstance 	};
enum vertex_cell_t			{ vertex_cell			};
enum vertex_dest_t			{ vertex_dest			};
enum vertex_idx_t			{ vertex_idx			};
enum vertex_mother_t		{ vertex_mother			};
enum vertex_deg_t			{ vertex_deg			};

enum edge_width_t			{ edge_width			};
enum edge_head_t			{ edge_head				};
enum edge_point_t			{ edge_point			};

namespace boost 
{
	BOOST_INSTALL_PROPERTY( vertex, diameter 		);
	BOOST_INSTALL_PROPERTY( vertex, width 			);
	BOOST_INSTALL_PROPERTY( vertex, coordinate 		);
	BOOST_INSTALL_PROPERTY( vertex, topoInstance 	);
	BOOST_INSTALL_PROPERTY( vertex, cell		 	);
	BOOST_INSTALL_PROPERTY( vertex, dest		 	);
	BOOST_INSTALL_PROPERTY( vertex, idx			 	);
	BOOST_INSTALL_PROPERTY( vertex, mother		 	);
	BOOST_INSTALL_PROPERTY( vertex, deg			 	);

	BOOST_INSTALL_PROPERTY( edge,   width		 	);
	BOOST_INSTALL_PROPERTY( edge,   head		 	);
	BOOST_INSTALL_PROPERTY( edge,   point		 	);
}
// boost geometry library deps
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/core/point_order.hpp> 
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/geometries/segment.hpp> 
// boost polygon library deps
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/voronoi_diagram.hpp>
#include <boost/polygon/voronoi_builder.hpp>
// fmt library deps
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#include "xlsxwriter.h"
#include <gnugraph.h>

// specialization for boost geometry to adapt to boost polygon point concept
#ifndef _BOOST_POLYGON_GEOMETRY_CONCEPT_POINT_XY
#define _BOOST_POLYGON_GEOMETRY_CONCEPT_POINT_XY
	namespace boost 
	{
		namespace polygon 
		{
			template <typename CoordinateType>
				struct geometry_concept<typename boost::geometry::model::d2::point_xy<CoordinateType>>
				{
					using type = typename boost::polygon::point_concept;
				};
		}
	}
#endif
#ifndef _BOOST_POLYGON_POINT_TRAITS_POINT_XY
#define _BOOST_POLYGON_POINT_TRAITS_POINT_XY
	namespace boost 
	{
		namespace polygon 
		{
			template <typename CoordinateType>
				struct point_traits<typename boost::geometry::model::d2::point_xy<CoordinateType>>
				{
					using coordinate_type = CoordinateType;
					static inline coordinate_type get( typename boost::geometry::model::d2::point_xy<CoordinateType> const& point, boost::polygon::orientation_2d orient )
					{
						return ( orient == HORIZONTAL ) ? boost::geometry::get<0>( point ) : boost::geometry::get<1>( point );
					}
				};
		}
	}
#endif
// register a vector of point as a geometry ring
BOOST_GEOMETRY_REGISTER_RING( std::vector<boost::polygon::point_data<double>> )
namespace boost::geometry
{
    template<>
    struct point_order<std::vector<boost::polygon::point_data<double>>>
    {
        static const order_selector value = counterclockwise;
        //static const order_selector value=clockwise;
    };
}

namespace syc::topology
{
	namespace detail::v_01							// implementation detail
	{
		// base_sketchable_forest
		template <typename coordinate_type, typename point_type>
			using vp = typename boost::property<boost::vertex_index_t,	int,												// reserved but not used
					   typename boost::property<vertex_diameter_t,		coordinate_type,									// reserved but not used
					   typename boost::property<vertex_width_t,			coordinate_type,									// reserved but not used
					   typename boost::property<vertex_coordinate_t, 	point_type,											// reserved but not used
					   typename boost::property<vertex_topoInstance_t,	std::array<std::list<std::vector<std::size_t>>, 2>,
					   typename boost::property<vertex_dest_t,			long int 
					   >>>>>>; 
		template <typename coordinate_type, typename point_type>
			using base_sketchable_forest = 
				typename boost::adjacency_list< 
					typename boost::vecS, 
					typename boost::vecS, 
					typename boost::bidirectionalS, 
					vp<coordinate_type, point_type>
				>;
		// topo_descr_impl
		using topo_descr_impl = 
			std::tuple<
				std::size_t, 									// vertex number this topoidentity belongs to
				int,											// 0: it belongs to an tree edge, 1: it belongs to an back edge, it belongs to a vertex 
				std::list<std::vector<std::size_t>>::iterator	// topoidentity positioin
			>;
		// slice_type
		class slice_type 
			: public std::list<topo_descr_impl>
		{
			// copy control
			public:
				slice_type() = default;
				slice_type( 
					bool 							d, 
					std::size_t 					s, 
					std::size_t 					t, 
					std::list<slice_type>::iterator p, 
					std::list<slice_type>::iterator n 
				): 
					direction( d ), 
					first( s ), 
					second( t ),
					pre( p ),
					nex( n ),
					mother( n )			// not robust
				{} 
			// data member
			public:
				std::size_t									index		=  1;		// reserved for algorithms, and reset_slice_indexes method. setting 1 is ungly!!
				bool										direction	=  0;		// 0: South, 1: North
				std::size_t									first		= -1; 
				std::size_t									second		= -1; 
				// for nets
				std::list<slice_type>::iterator				pre;
				std::list<slice_type>::iterator				nex;
				// for slice hierarchy (with no default initializations)
				std::list<slice_type>::iterator				mother;
				std::list<std::list<slice_type>::iterator>	children;
		};
		// initialize_topology	 
		template <typename sketable_forest>
		class initialize_topology : 
			public boost::default_dfs_visitor 
		{
			public:
				initialize_topology( sketable_forest& f ): _F( f ) {}

				template <typename Edge, typename Graph>	
					void examine_edge( Edge e, Graph const& g )  
					{
						auto u = boost::source( e, g );
						auto& Tu = boost::get( vertex_topoInstance, _F, u );
						auto p = Tu[ 1 ].insert( Tu[ 1 ].end(), std::vector<std::size_t>( 1, _topocnt++ ) );
						_F.sketch.front().emplace_back( u, 2, p );
						auto v = boost::target( e, g );
						auto& Tv = boost::get( vertex_topoInstance, _F, v );
						auto q = Tv[ 0 ].insert( Tv[ 0 ].end(), std::vector<std::size_t>( 1, _topocnt++ ) );
						_F.sketch.front().emplace_back( v, 0, q );
					}
				template <typename Vertex, typename Graph>	
					void finish_vertex( Vertex u, Graph const& g )  
					{
						auto& T = boost::get( vertex_topoInstance, _F, u );
						if ( boost::in_degree( u, _F ) == 0 ) { 	
							auto p = T[ 1 ].insert( T[ 1 ].end(), std::vector<std::size_t>( 1, _topocnt++ ) );
							_F.sketch.front().emplace_back( u, 2, p );
						}
						else { 
							auto p = T[ 1 ].insert( T[ 1 ].end(), std::vector<std::size_t>( 1, _topocnt++ ) );
							_F.sketch.front().emplace_back( u, 2, p );
							T[ 0 ].front().push_back( _topocnt++ );
							_F.sketch.front().emplace_back( u, 1, T[ 0 ].begin() );
						}
					}
			private:
				sketable_forest&		 	_F; 
				std::size_t					_topocnt = 0;
		};
		// compute_dest 
		template <typename sketable_forest>
		class compute_dest: public boost::default_bfs_visitor 
		{
			// copy control
			public:
				compute_dest( sketable_forest& f, unsigned i, unsigned n ): _F( f ), _i( i ), _num_roots( n ) {}
			// data member
			private:
				sketable_forest& _F; 
				unsigned _i;
				long int _num_roots;
			// events
			public:	
				// tree_edge
				template <typename Edge, typename Graph>	
				void tree_edge( Edge e, Graph const& g )  
				{
					auto u = boost::source( e, g ); 
					auto v = boost::target( e, g ); 

					auto du = u;
					if ( 0 != boost::in_degree( du, _F ) ) {				// it is not a root
						du = *boost::adjacent_vertices( du, _F.nets ).first;
					}
					auto dv = v;
					if ( 0 != boost::in_degree( dv, _F ) ) {				// it is not a root
						dv = *boost::adjacent_vertices( v, _F.nets ).first;
					}
					long int dis = dv - du;
					long int dis_prime = ( dis <= 0 )? _num_roots + dis: -_num_roots + dis; 
					auto dest_u = boost::get( vertex_dest, _F, u );
					dis += dest_u;
					dis_prime += dest_u;
					auto x = ( abs( dis ) < abs( dis_prime ) )? dis: dis_prime;
					boost::put( vertex_dest, _F, v, x ); 
				}
		};
		// net_type_impl
		template <typename coordinate_type>
		using net_type_impl = 
			boost::adjacency_list<
				boost::vecS, 
				boost::vecS, 
				boost::bidirectionalS, 
				boost::no_property, 
				boost::property<boost::edge_index_t, std::size_t,						// reserved but not used
					boost::property<boost::edge_name_t, std::string,					// reserved but not used
						boost::property<edge_width_t, coordinate_type,					// reserved but not used
							boost::property<edge_head_t, typename slice_type::iterator	
				>>>>
			>;
		// delaunay_triangulation_graph_type
		template <typename CoordinateType>
		using delaunay_triangulation_graph_type = 
			boost::adjacency_list<
				boost::vecS, 
				boost::vecS, 
				boost::undirectedS,
				boost::property<boost::vertex_index_t, std::size_t,
					boost::property<vertex_cell_t, typename boost::polygon::voronoi_diagram<CoordinateType>::const_cell_iterator
				>>,
				boost::property<boost::edge_weight_t, CoordinateType,
					boost::property<boost::edge_index_t, std::size_t,			// Used as a marker. 0: Triangulation Graph Edge 1: MST Edge
						boost::property<boost::edge_color_t, int,				// Used as a tempary marker for traversing
							boost::property<boost::edge_capacity_t, unsigned 	// Used as a recording of # of traces passing an edge
				>>>>
			>;
		// build_sketchable_forest
		template <typename delaunay_triangulation_graph_type, typename sketable_forest>
		class build_sketchable_forest: public boost::default_bfs_visitor 
		{
			// copy control
			public:
				build_sketchable_forest( delaunay_triangulation_graph_type& dgt, sketable_forest& sf ): _DTG( dgt ), _SF( sf ) {}
			// data member
			private:
				delaunay_triangulation_graph_type&	_DTG;
				sketable_forest&					_SF;
			// events
			public:	
				// examine_vertex
				template <typename Vertex, typename Graph>	
					void examine_vertex( Vertex u, Graph const& g )  
					{
						// reset color of all adjacent edges
						for ( auto ei = boost::out_edges( u, _DTG ); ei.first !=  ei.second; ++ei.first ) {
							boost::put( boost::edge_color, _DTG, *ei.first, 0 );
						}
					}
				// tree_edge
				template <typename Edge, typename Graph>	
					void tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, _DTG, e, 1 );
					}
				// non_tree_edge
				template <typename Edge, typename Graph>	
					void non_tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, _DTG, e, 2 );
					}
				// finish_vertex
				template <typename Vertex, typename Graph>	
					void finish_vertex( Vertex u, Graph const& g )  
					{
						auto it_end	= boost::get( vertex_cell, _DTG, u )->incident_edge();
						auto it		= it_end->prev();
						while ( it != it_end ) {
							auto v = it->twin()->cell()->source_index();
							auto e = boost::edge( u, v, _DTG ).first;
							if ( 2 == boost::get( boost::edge_color, _DTG, e ) ) break;
							it = it->prev();
						}
						auto it2 = it;
						do {
							auto v = it2->twin()->cell()->source_index();
							auto e = boost::edge( u, v, _DTG ).first;
							if ( 1 == boost::get( boost::edge_color, _DTG, e ) ) {
								boost::add_edge( u, v, _SF );
							}
							it2 = it2->prev();
						} while ( it2 != it );
					}
		};
		// is_MST_edge
		template <typename EdgeIndexMap>
		struct is_MST_edge 
		{
			public:
				is_MST_edge() = default;
				is_MST_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 1;
				}
			private:
				EdgeIndexMap I;
		};
	}
	namespace model::v_01								// model for users
	{
		// sketchable_forest
		template <typename coordinate_type = double, typename point_type = boost::polygon::point_data<coordinate_type>>
			class sketchable_forest
				: public syc::topology::detail::v_01::base_sketchable_forest<coordinate_type, point_type>
			{
				// friend declaration
				template <typename T> 
					friend class syc::topology::detail::v_01::initialize_topology;
				// type declaration
				public:
					using net_type				= detail::v_01::net_type_impl<coordinate_type>;
					using slice_type 			= detail::v_01::slice_type;
					using sketch_type 			= std::list<slice_type>;
					using topo_descriptor 		= typename slice_type::iterator;
					using slice_descriptor		= typename sketch_type::iterator;
					using topology_type			= std::list<std::pair<slice_descriptor,topo_descriptor>>;
					using result_type			= std::vector<std::list<std::pair<std::size_t, int>>>;
					enum _SliceDir	{ downward = 0, upward = 1 };
					enum _Att		{ tree_side_of_an_edge = 0, back_side_of_an_edge = 1, vertex = 2 };
				// member function
				public:
					// initialization operation
					void											make_sketch( void ); 			
					void											reset_slice_indexes( void );
					// debug operations
					void											pretty_printer( void ) const;
					void											txt_vitualizer( void ) const; 
					void											txt_vitualizer2( std::ostream& os ) const; 
					void											print_all_2_pin_nets( void ) const; 
					// proviided operations for vertex descriptor v
					std::vector<topo_descriptor>					topo_instances( std::size_t const v );
					// provided for topo_descrpitor
					std::size_t										corresponding_vertex( topo_descriptor const d );
					int												attribution( topo_descriptor const d );
					std::size_t										topo_id( topo_descriptor const d );								
					topo_descriptor									ordering_pair( topo_descriptor const d );
					slice_descriptor								slice( topo_descriptor const d );	
                    bool                                            shortest_direction( topo_descriptor s, topo_descriptor t );
					topo_descriptor									make_slice( topo_descriptor const s, topo_descriptor const t );
					topo_descriptor									free_slice( topo_descriptor const arrow );
					result_type										results( void );
					// provided operations for topo-instance s, t
					std::size_t										make_slice( std::size_t const s, std::size_t const t );
					std::size_t										free_slice( std::size_t const arrow );
				private:
					// provided operations for results
					void											make_slice_hierarchy( void );
				// data member
				public:
					net_type		nets;						// reserved only for two pin nets
				private:
					sketch_type 	sketch;
					topology_type 	topology;

					std::vector<typename topology_type::iterator>	
						position;

					std::vector<
						std::pair<
							slice_descriptor, 
							bool								// 0: start, 1: end
					>>	attachment;
				// coppy control
				public:
					sketchable_forest(): 
						sketch( 1 ) 
					{}
					sketchable_forest( std::size_t n ): 
						syc::topology::detail::v_01::base_sketchable_forest<coordinate_type, point_type>( n ),
						sketch( 1 ) 
					{}
			};
		// sketchable_forest make_sketch
		template <typename coordinate_type, typename point_type>
			inline void sketchable_forest<coordinate_type, point_type>::
			make_sketch( void ) 
			{ 
				boost::depth_first_search( *this, visitor( syc::topology::detail::v_01::initialize_topology<decltype(*this)>(*this) ) ); 
				auto it = sketch.begin(), it_end = sketch.end();
				auto it2 = it->begin();
				std::size_t n = it->size();
				for ( std::size_t i = 0; i != n; ++i, ++it2 ) {
					auto p = topology.emplace( topology.end(), it, it2 );
					position.push_back( p );
				}
				attachment.resize( n, std::pair<sketch_type::iterator, bool>( it_end, 0) );
				sketch.front().pre = it_end;
				sketch.front().nex = it_end;
				sketch.front().mother = it_end;
			}
		// sketchable_forest reset_slice_indexes
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			reset_slice_indexes( void ) -> void
			{
				for( auto& slice: sketch ) slice.index = 0; 
			}
		// sketchable_forest pretty_printer
		template <typename coordinate_type, typename point_type>
			inline void sketchable_forest<coordinate_type, point_type>::
			pretty_printer( void ) const 
			{
                /*
				for ( std::size_t i = 0, n = boost::num_vertices( *this ); i != n; ++i ) {
					std::cout << "vertex: " << i << std::endl;

					std::cout << "\tForest Edge List: ";
					for ( auto vi = boost::adjacent_vertices( i, *this ); vi.first != vi.second; ++vi.first ) {
						std::cout << *vi.first << " ";
					}
					std::cout << std::endl;

					std::cout << "\tNets Edge List: ";
					for ( auto vi = boost::adjacent_vertices( i, nets ); vi.first != vi.second; ++vi.first ) {
						std::cout << *vi.first << " ";
					}
					std::cout << std::endl;

					std::cout << "\tdiameter: " 	<< boost::get( vertex_diameter, *this, i ) 								<< std::endl;
					std::cout << "\twidth: " 		<< boost::get( vertex_width, *this, i ) 								<< std::endl;
					std::cout << "\tcoordinate: " 	<< boost::geometry::wkt( boost::get( vertex_coordinate, *this, i ) ) 	<< std::endl;

					auto T = boost::get( vertex_topoInstance, *this, i );
					std::cout << "\ttopoInstances[ 0 ]: ";
					for ( auto const& x: T[ 0 ] ) {
						std::cout << "(" << x[0] << "," << x[1] << ") ";
					}
					std::cout << std::endl;
					std::cout << "\ttopoInstances[ 1 ]: ";
					for ( auto const& x: T[ 1 ] ) {
						std::cout << x[0] << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "---" << std::endl;
				std::string A = "TopoInstance", B = " position", C = " TopoF", D = " TopoS", E = " attachF", F = " attachS";
				std::cout << std::setw( A.size() ) << A
						  << std::setw( B.size() ) << B
						  << std::setw( C.size() ) << C
						  << std::setw( D.size() ) << D
						  << std::setw( E.size() ) << E
						  << std::left
						  << std::setw( F.size() ) << F 
						  << std::endl;
				std::cout << "=======================================================================" << std::endl;
				for ( std::size_t i = 0, i_end = position.size(); i != i_end; ++i ) {
					std::cout << std::setw( A.size() ) << i
//							  << std::setw( B.size() ) << std::distance( topology.begin(), position[ i ] ) 						
//							  << std::setw( C.size() ) << std::distance( sketch.begin(), position[ i ]->first )					
//							  << std::setw( D.size() ) << std::distance( position[ i ]->first->begin(), position[ i ]->second ) 
//							  << std::setw( E.size() ) << std::distance( sketch.begin(), attachment[ i ].first )				
							  << std::boolalpha
							  << std::setw( F.size() ) << attachment[ i ].second
							  << std::endl;
				}
				std::cout << "---" << std::endl;
                */
				std::cout << "sketch" << std::endl;
				std::cout << "=======================================================================" << std::endl;
				std::size_t i = 0;
				for ( auto const& sl: sketch ) {
					std::cout << "slice: " 	<< i++ 			<< std::endl;
					std::cout << "\tmembers: ";
					for ( auto const& td: sl ) {
						std::cout << "(" << std::get<0>( td ) << "," << std::get<1>( td ) << ",[";
						for ( auto const& e: *std::get<2>( td ) ) {
							std::cout << " " << e; 
						}
						std::cout << " ]) ";
					}
					std::cout << std::endl;
					std::cout << "\tindex:     " 	<< sl.index 		<< std::endl;
					std::cout << "\tdirection: " 	<< sl.direction 	<< std::endl;
					std::cout << "\tfirst:     " 	<< sl.first 		<< std::endl;
					std::cout << "\tsecond:    " 	<< sl.second 		<< std::endl;
//					std::cout << "\tpre:       " 	<< std::distance( sketch.begin(), sl.pre ) << std::endl;
//					std::cout << "\tnex:       " 	<< std::distance( sketch.begin(), sl.nex ) << std::endl;
//					std::cout << "\tmother:    " 	<< std::distance( sketch.begin(), sl.mother ) << std::endl;
					std::cout << "\tchildren:  ";
					for ( auto const& i : sl.children ) 
                    {
//						 std::cout << std::distance( sketch.begin(), i ) << " ";
                    }
					std::cout << std::endl;

				}
			}
		// sketchable_forest txt_vitualizer
		template <typename coordinate_type, typename point_type>
			inline void sketchable_forest<coordinate_type, point_type>::
			txt_vitualizer( void ) const  
			{
				// prepare
				unsigned 					n = sketch.size() - 1, 
											m = topology.size();
				unsigned					num_width = log10(m) + 2;
				std::vector<std::string> 	R( m, std::string( n, ' ' ) );
				std::size_t 				j = 0;
				std::vector<std::size_t> 	index( m );
				for ( auto const p: topology ) {
					auto& t = *p.second;
					std::size_t i = ( std::get<1>( t ) == 1 ) ? (*std::get<2>( t ))[ 1 ] : (*std::get<2>( t ))[ 0 ];
					index[ i ] = j++;
				}
				// initialize R
				j = 0;
				for ( auto it = sketch.begin(), it_end = sketch.end(); ++it != it_end; ++j ) {
					std::size_t s = index[ it->first ];
					std::size_t t = index[ it->second ];
					std::size_t k1 = 0, k2 = 0;
					while ( k1 <= j ) 
						R[ s ][ k1++ ] = '-';
					while ( k2 <= j )
						R[ t ][ k2++ ] = '-';
					R[ s ][ 0 ] = '>';
					R[ t ][ 0 ] = '<';
					if ( it->direction == 0 ) {
						R[ s ][ j ] = '+';
						while ( ++s < t ) {
							R[ s ][ j ] = ( it->direction ) ? '|' : '!';
						}
						R[ s ][ j ] = '+';
					}
					else {
						R[ s ][ j ] = '+';
						while ( --s > t ) {
							R[ s ][ j ] = ( it->direction ) ? '|' : '!';
						}
						R[ s ][ j ] = '+';
					}
				}
				// output
				auto it = topology.begin();
				std::cout << std::string( num_width * 2 + 4 + n + 1, '=' ) << std::endl;
				for ( auto const& row: R ) {
					auto& t = *it->second;
					std::size_t i = ( std::get<1>( t ) == 1 ) ? (*std::get<2>( t ))[ 1 ] : (*std::get<2>( t ))[ 0 ];
					if 		( std::get<1>( t ) == 0 ) {
						std::cout << std::setw( num_width ) << " "
								  << " | +" 
								  << std::right
								  << std::setw( num_width ) << i << " ";
					}
					else if ( std::get<1>( t ) == 1 ) {
						std::cout << std::setw( num_width ) << " "
								  << " | -" 
								  << std::right
								  << std::setw( num_width ) << i << " "; 
					}
					else { 
						std::cout << std::setw( num_width ) << std::get<0>( t )
								  << " |  "
								  << std::right
								  << std::setw( num_width ) << i << " "; 
					}
					std::cout << row << std::endl;
					++it;
				}
				std::cout << std::string( num_width * 2 + 4 + n + 1, '=' ) << std::endl;
			}
		// sketchable_forest txt_vitualizer2
		template <typename coordinate_type, typename point_type>
			inline void sketchable_forest<coordinate_type, point_type>::
			txt_vitualizer2( std::ostream& os ) const  
			{
				// prepare
				unsigned 					n = sketch.size() - 1, 
											m = topology.size();
				unsigned					num_width = log10(m) + 2;
				std::vector<std::string> 	R( m, std::string( n, ' ' ) );
				std::size_t 				j = 0;
				std::vector<std::size_t> 	index( m );
				for ( auto const p: topology ) {
					auto& t = *p.second;
					std::size_t i = ( std::get<1>( t ) == 1 ) ? (*std::get<2>( t ))[ 1 ] : (*std::get<2>( t ))[ 0 ];
					index[ i ] = j++;
				}
				// initialize R
				j = 0;
				for ( auto it = sketch.begin(), it_end = sketch.end(); ++it != it_end; ++j ) {
					std::size_t s = index[ it->first ];
					std::size_t t = index[ it->second ];
					std::size_t k1 = 0, k2 = 0;
					while ( k1 <= j ) 
						R[ s ][ k1++ ] = '-';
					while ( k2 <= j )
						R[ t ][ k2++ ] = '-';
					R[ s ][ 0 ] = '>';
					R[ t ][ 0 ] = '<';
					if ( it->direction == 0 ) {
						R[ s ][ j ] = '+';
						while ( ++s < t ) {
							R[ s ][ j ] = ( it->direction ) ? '|' : '!';
						}
						R[ s ][ j ] = '+';
					}
					else {
						R[ s ][ j ] = '+';
						while ( --s > t ) {
							R[ s ][ j ] = ( it->direction ) ? '|' : '!';
						}
						R[ s ][ j ] = '+';
					}
				}
				// output
				auto it = topology.begin();
				os << std::string( num_width * 2 + 4 + n + 1, '=' ) << std::endl;
				for ( auto const& row: R ) {
					auto& t = *it->second;
					std::size_t i = ( std::get<1>( t ) == 1 ) ? (*std::get<2>( t ))[ 1 ] : (*std::get<2>( t ))[ 0 ];
					if 		( std::get<1>( t ) == 0 ) {
						os << std::setw( num_width ) << " "
								  << " : +" 
								  << std::right
								  << std::setw( num_width ) << i << " ";
					}
					else if ( std::get<1>( t ) == 1 ) {
						os << std::setw( num_width ) << " "
								  << " : -" 
								  << std::right
								  << std::setw( num_width ) << i << " "; 
					}
					else { 
						os << std::setw( num_width ) << std::get<0>( t )
								  << " :  "
								  << std::right
								  << std::setw( num_width ) << i << " "; 
					}
					os << row << std::endl;
					++it;
				}
				os << std::string( num_width * 2 + 4 + n + 1, '=' ) << std::endl;
			}
		// sketchable_forest print_all_2_pin_nets
		template <typename coordinate_type, typename point_type>
		inline void sketchable_forest<coordinate_type, point_type>::
		print_all_2_pin_nets( void ) const 
		{
			unsigned i = 0;
			for ( auto ei = boost::edges( nets ); ei.first != ei.second; ++ei.first, ++i )
			{
				std::cout << "net " << i << ": (" << boost::source( *ei.first, nets ) << ", " << boost::target( *ei.first, nets ) << ")\n";
			}
		}
		// sketchable_forest topo_instances
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			topo_instances( std::size_t const v ) -> std::vector<topo_descriptor>
			{ 
				std::vector<topo_descriptor> R;
				auto const& T = boost::get( vertex_topoInstance, *this, v );
				for( auto const& x: T[ 1 ] ) {
					R.emplace_back( position[ x[ 0 ] ]->second );
				}
				return R;
			}
		// sketchable_forest corresponding_vertex
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			corresponding_vertex( topo_descriptor const d )	-> std::size_t
			{ 
				return std::get<0>( *d ); 
			}
		// sketchable_forest attribution
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			attribution( topo_descriptor const d ) -> int							
			{ 
				return std::get<1>( *d ); 
			}
		// sketchable_forest topo_id
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			topo_id( topo_descriptor const d ) -> std::size_t								
			{ 
				return ( std::get<1>( *d ) == 1 ) ? (*std::get<2>( *d ))[ 1 ] : (*std::get<2>( *d ))[ 0 ]; 
			}
		// sketchable_forest ordering_pair
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			ordering_pair( topo_descriptor const d ) -> topo_descriptor
			{ 
				auto attr = std::get<1>( *d ); 
				if		( attr == 0 ) {
					return position[ (*std::get<2>( *d ))[ 1 ] ]->second;
				}
				else if ( attr == 1 ) {
					return position[ (*std::get<2>( *d ))[ 0 ] ]->second;
				}
				else {
					std::cout << "error: cannot find ordering pair from a topodescriptor of a vertex atrribution" << std::endl;
				}
			}
		// sketchable_forest slice
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			slice( topo_descriptor const d ) -> slice_descriptor	
			{ 
				return ( ( std::get<1>( *d ) == 1 ) ? position[ (*std::get<2>( *d ))[ 1 ] ]->first : position[ (*std::get<2>( *d ))[ 0 ] ]->first );
			}
		// sketchable_forest shortest_direction
		template <typename coordinate_type, typename point_type>
			bool sketchable_forest<coordinate_type, point_type>::
			shortest_direction( topo_descriptor s, topo_descriptor t )  // CW: false. CCW: true.
			{ 
		        auto ori_mark = this->slice( s )->end();
                auto ccw_it = s;
                auto cw_it  = s;
                while ( true )
                {
                    if ( ++ccw_it == ori_mark ) ++ccw_it;
                    if ( --cw_it  == ori_mark ) --ccw_it;
                    if ( ccw_it == t ) return true;
                    if ( cw_it  == t ) return false;
                }
			}
		// sketchable_forest make_slice 
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			make_slice( topo_descriptor const s, topo_descriptor const t ) -> topo_descriptor
			{
				std::size_t i, j;
				auto attr_s = std::get<1>( *s ); 
				auto attr_t = std::get<1>( *t ); 
				if		( attr_s == 0 ) {
					i = (*std::prev(std::get<2>( *s )))[ 0 ];
				}
				else if ( attr_s == 1 ) {
					i = (*std::prev(std::get<2>( *s )))[ 1 ];
				}
				else {
					i = (*std::get<2>( *s ))[ 0 ];
				}
				j = ( attr_t == 1 ) ? (*std::get<2>( *t ))[ 1 ] : (*std::get<2>( *t ))[ 0 ];
				return position[ make_slice( i, j ) ]->second;
			}
		// sketchable_forest free_slice 
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			free_slice( topo_descriptor const arrow ) -> topo_descriptor
			{
				std::size_t i;
				auto attr = std::get<1>( *arrow ); 
				if		( attr == 0 ) {
					i = (*std::prev(std::get<2>( *arrow )))[ 0 ];
				}
				else if ( attr == 1 ) {
					i = (*std::prev(std::get<2>( *arrow )))[ 1 ];
				}
				else {
					i = (*std::get<2>( *arrow ))[ 0 ];
				}
				return position[ free_slice( i ) ]->second;
			}
		// sketchable_forest make_slice
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			make_slice( std::size_t const s, std::size_t const t ) -> std::size_t 
			{
				// copy info of the source and target
				std::size_t	n 			= position.size();

				auto 		pos_t 		= position[ t ];
				auto 		sldescr_t	= pos_t->first;
				auto 		tdescr_t	= pos_t->second;
				auto 		vertex_t	= std::get<0>( *tdescr_t );
				auto 		attr_t 		= std::get<1>( *tdescr_t );
				auto 		tinsts_t 	= std::get<2>( *tdescr_t );

				auto 		pos_s 		= position[ s ];
				auto 		sldescr_s	= pos_s->first;
				auto 		tdescr_s	= pos_s->second;
				auto 		vertex_s	= std::get<0>( *tdescr_s );
				auto 		attr_s 		= std::get<1>( *tdescr_s );
				auto 		tinsts_s 	= std::get<2>( *tdescr_s );

				bool 		direction	= 0;		
					while ( pos_s != pos_t )  {
						if 		( pos_s == topology.end() ) {
							direction = 1;
							break;
						}
						++pos_s;
					}; 
					pos_s = position[ s ];

				auto 		h 			= sketch.insert( std::next( sldescr_t ), slice_type( direction, s, t, attachment[ s ].first, sketch.end() ) );

				auto 		c1 			= ( direction == 0 ) ? tdescr_s : tdescr_t;
				auto 		c2 			= ( direction == 0 ) ? tdescr_t : tdescr_s;

				int  		i			= 0;
				int  		j			= 0;
				// operations
				auto add_edge_for_target = [&]() -> void {
					auto& 		T = boost::get( vertex_topoInstance, *this, vertex_t );
					std::vector<std::size_t> temp( 2, n++ );
					h->second = temp[ i ];
					++temp[ i ];
					auto p = T[ 0 ].insert( tinsts_t, std::move( temp ) );
					auto q = sldescr_t->emplace( tdescr_t, vertex_t, attr_t, p ); 
					auto m = topology.emplace( pos_t, sldescr_t, q );
					position.push_back( m );
					attachment.emplace_back( h, 1 );
					position.push_back( position[ (*tinsts_t)[ i ] ] );
					attachment.emplace_back( h, 0 );
				};
				auto add_edge_for_source = [&]() -> void {
					auto q = sldescr_s->emplace( tdescr_s, vertex_s, attr_s, tinsts_s ); 
					auto m = topology.emplace( pos_s, sldescr_s, q );
					position[ s ] = m;
					attachment[ s ].first = h;
					attachment[(*tinsts_s)[j]].first->nex = h;
				};
				auto splice_to_new_slice = [&]() -> void {
					auto prev_h = std::prev( h );
					h->splice( h->end(), *prev_h, prev_h->begin(), c1            );
					h->splice( h->end(), *prev_h, c2,              prev_h->end() );
					for ( auto it = h->begin(), it_end = h->end(); it != it_end; ++it ) {
						std::size_t inst = ( std::get<1>( *it ) == 1 )? (*std::get<2>( *it ))[ 1 ]: (*std::get<2>( *it ))[ 0 ];
						position[ inst ]->first = h;
						position[ inst ]->second = it;
					}
				};


				// check target attribution
				if 		( attr_t == 0 ) { 		// case 1 and case 2
					i = 1;
					add_edge_for_target();
				}
				else if ( attr_t == 1 ) {		// case 3 and case 4
					++tdescr_t;			
					++pos_t;
					i = 0;
					add_edge_for_target();
					if ( direction == 0 ) 	++c2;
					else					++c1;
				}
				// check source attribution
				if 		( attr_s == 1 ) {		// case 1 and case 2	
					j = 0;
					--tinsts_s;
					++tdescr_s;
					++pos_s;
					add_edge_for_source();
					if ( direction == 0 ) 	++c1;
					else					++c2;
				}
				else if ( attr_s == 0 ) { 		// case 3 and case 4
					j = 1;
					--tinsts_s;
					add_edge_for_source();
				}

				if ( attr_t == 2 ) { 		
					attachment[ t ].first 	= h;
					attachment[ t ].second 	= 1;
					if ( direction == 1 ) 	++c1;
					n = t;							// for return value
				}
				if ( attr_s == 2 ) { 		
					attachment[ s ].first 	= h;
					if ( direction == 0 ) 	++c1;
				}

				// splice
				splice_to_new_slice();
				return n;
			}
		// sketchable_forest free_slice
		template <typename coordinate_type, typename point_type>
			inline auto sketchable_forest<coordinate_type, point_type>::
			free_slice( std::size_t const arrow ) -> std::size_t 
			{
				if ( attachment[ arrow ].first == sketch.end() ) return arrow;
				// info
				auto new_slice	= attachment[ arrow ].first;
				auto old_slice	= std::prev( new_slice );

				auto s 			= new_slice->first;
				auto pos_s 		= position[ s ];
				auto sldescr_s	= pos_s->first;
				auto tdescr_s	= pos_s->second;
				auto vertex_s	= std::get<0>( *tdescr_s );
				auto attr_s 	= std::get<1>( *tdescr_s );
				auto tinsts_s 	= std::get<2>( *tdescr_s );

				auto t 			= new_slice->second;
				auto pos_t 		= position[ t ];
				auto sldescr_t	= pos_t->first;
				auto tdescr_t	= pos_t->second;
				auto vertex_t	= std::get<0>( *tdescr_t );
				auto attr_t 	= std::get<1>( *tdescr_t );
				auto tinsts_t 	= std::get<2>( *tdescr_t );

				auto pos_c2 	= pos_t; 
				std::size_t i 	= 0;
				// operation
				auto splice_to_prev_slice = [&]() -> void {
					auto it_end = old_slice->begin();
					auto it2	= std::prev( old_slice->end() );
					old_slice->splice( old_slice->begin(), *new_slice, new_slice->begin(), pos_c2->second   );
					old_slice->splice( old_slice->end(),   *new_slice, pos_c2->second,     new_slice->end() );
					for ( auto it = old_slice->begin(); it != it_end; ++it ) {
						std::size_t inst = ( std::get<1>( *it ) == 1 )? (*std::get<2>( *it ))[ 1 ]: (*std::get<2>( *it ))[ 0 ];
						position[ inst ]->first = old_slice;
						position[ inst ]->second = it;
					}
					for ( auto it2_end = old_slice->end(); ++it2 != it2_end; ) {
						std::size_t inst = ( std::get<1>( *it2 ) == 1 )? (*std::get<2>( *it2 ))[ 1 ]: (*std::get<2>( *it2 ))[ 0 ];
						position[ inst ]->first = old_slice;
						position[ inst ]->second = it2;
					}
					if ( new_slice->pre != sketch.end() ) new_slice->pre->nex = sketch.end();
				};
				auto erase_edge_for_target = [&]() -> void {
					auto&	T = boost::get( vertex_topoInstance, *this, vertex_t );
					T[ 0 ].erase( tinsts_t );
					sldescr_t->erase( tdescr_t ); 
					topology.erase( pos_t );
					position.pop_back();
					position.pop_back();
					attachment.pop_back();
					attachment.pop_back();
				};
				auto erase_edge_for_source = [&]() -> void {
					attachment[ s ].first = attachment[ s ].first->pre;
					std::advance( position[ s ], i );
					sldescr_s->erase( tdescr_s ); 
					topology.erase( pos_s );
				};
				// check direction 
				if ( new_slice->direction == 0 ) {		
					// check target attribution
					if 		( attr_t == 0 ) {
						++pos_c2;
						splice_to_prev_slice();
						erase_edge_for_target();
					}
					else if ( attr_t == 1 ) {
						splice_to_prev_slice();
						erase_edge_for_target();
					}
					else {
						splice_to_prev_slice();
						attachment[ t ].first = sketch.end();
						attachment[ t ].second = 0;
					}	
					// check source attribution
					if 		( attr_s == 0 ) {
						++i;
						erase_edge_for_source();
					}
					else if ( attr_s == 1 ) {
						--i;
						erase_edge_for_source();
					}
					else {
						attachment[ s ].first = sketch.end();
					}	
					sketch.erase( new_slice );
				}
				else {
					pos_c2 = position[ s ];
					// check source attribution
					if 		( attr_s == 0 ) {
						++pos_c2;
						splice_to_prev_slice();
						++i;
						erase_edge_for_source();
					}
					else if ( attr_s == 1 ) {
						--i;
						splice_to_prev_slice();
						erase_edge_for_source();
					}
					else {
						splice_to_prev_slice();
						attachment[ s ].first = sketch.end();
					}	
					// check target attribution
					if 		( attr_t == 0 ) {
						erase_edge_for_target();
					}
					else if ( attr_t == 1 ) {
						erase_edge_for_target();
					}
					else {
						attachment[ t ].first = sketch.end();
						attachment[ t ].second = 0;
					}	
					sketch.erase( new_slice );
				}
				return s;
			}
		// sketchable_forest make_slice_hierarchy
		template <typename coordinate_type, typename point_type>
		inline void sketchable_forest<coordinate_type, point_type>::
		make_slice_hierarchy( void ) 
		{
			unsigned i = 0;
			std::vector< std::size_t > frame_index;		// map topo_id to it's index on the frame
			for ( unsigned i = 0; i != position.size(); ++i ) {
				frame_index.push_back( std::distance( topology.begin(), position[ i ] ) );
			}
			auto South 	= [&]( slice_type const& s1, slice_type const& s2 )-> bool {
				return ( std::min( frame_index[s1.first], frame_index[s1.second] ) > std::max( frame_index[s2.first], frame_index[s2.second] ) );
			}; 
			auto North 	= [&]( slice_type const& s1, slice_type const& s2 )-> bool {
				return ( std::max( frame_index[s1.first], frame_index[s1.second] ) < std::min( frame_index[s2.first], frame_index[s2.second] ) );
			}; 
			auto East 	= [&]( slice_type const& s1, slice_type const& s2 )-> bool {
				return ( std::min( frame_index[s1.first], frame_index[s1.second] ) < std::min( frame_index[s2.first], frame_index[s2.second] ) && 
						 std::max( frame_index[s1.first], frame_index[s1.second] ) > std::max( frame_index[s2.first], frame_index[s2.second] )    );
			}; 
			auto root = sketch.begin();
			root->mother = root;
			for ( auto root = sketch.begin(), slice_it = root; ++slice_it != sketch.end(); ) {
				auto it = root->children.begin(), it_end = root->children.end();
				for ( ;it != it_end; ++it ) {
					if ( North( *slice_it, **it ) ) break;
					if ( South( *slice_it, **it ) ) continue;
					auto m_it = it;
					while ( it != it_end && East( *slice_it, **it )  ) {
						(*it)->mother = slice_it;
						++it;
					}
					slice_it->children.splice( slice_it->children.end(), root->children, m_it, it ); 
					break;
				}
				root->children.insert( it, slice_it );
				slice_it->mother = root;
			}
		}
		// sketchable_forest results( void )
		template <typename coordinate_type, typename point_type>
		inline auto sketchable_forest<coordinate_type, point_type>::
		results( void ) -> std::vector<std::list<std::pair<std::size_t, int>>>
		{
			std::vector<std::list<std::pair<std::size_t, int>>> r;
			make_slice_hierarchy();			
			for ( auto ni = boost::edges( nets ); ni.first != ni.second; ++ni.first ) 
			{
				std::list<std::pair<std::size_t, int>> path;
				// finding s as the only topo_descriptor attached to a slice_descriptor and being the beginning of the net.
				auto vertex_s = boost::source( *ni.first, nets ); 
				auto const& T = boost::get( vertex_topoInstance, *this, vertex_s );
				auto s = position[ T[1].front()[0] ]->second;
				for ( auto it = T[1].begin(); ++it != T[1].end(); ) 
				{
					if ( attachment[ (*it)[0] ].first != sketch.end() )
					{
						s = position[ (*it)[0] ]->second;
					}
				}
				for ( auto slice_it = attachment[ topo_id( s ) ].first; slice_it != sketch.end(); slice_it = slice_it->nex ) {
					auto a1 = position[ slice_it->first  ];
					auto a2 = position[ slice_it->second ];
					if ( slice_it->direction == 0 ) { 
						for ( ; a1 != a2; ++a1 ) {
							if ( attribution( a1->second ) == 2 ) {		
								auto v = corresponding_vertex( a1->second );
								auto x = attachment[ topo_id( a1->second ) ].first;
								int i = 0;
								if ( x == sketch.end() ) {
									x = std::next( a1->first );
									++i;
								}
								while ( x != slice_it ) {
									++i;
									x = x->mother;
								}
								path.emplace_back( v, i );
							}
						}
						if ( attribution( a2->second ) == 2 ) {		
							auto v = corresponding_vertex( a2->second );
							path.emplace_back( v, 0 );
						}
					}
					else {
						for ( ; a1 != a2; --a1 ) {
							if ( attribution( a1->second ) == 2 ) {		
								auto v = corresponding_vertex( a1->second );
								auto x = attachment[ topo_id( a1->second ) ].first;
								int i = 0;
								if ( x == sketch.end() ) {
									x = std::next( a1->first );
									--i;
								}
								while ( x != slice_it ) {
									--i;
									x = x->mother;
								}
								path.emplace_back( v, i );
							}
						}
						if ( attribution( a2->second ) == 2 ) {		
							auto v = corresponding_vertex( a2->second );
							path.emplace_back( v, 0 );
						}
					}
				}
				if ( !path.empty() ) // then path contains at least 2 elements
				{ 
					for ( auto it1 = path.begin(), it2 = std::next( it1 ); it2 != path.end(); ++it2 )
					{
						if ( it1->first == it2->first ) 
							it1 = path.erase( it1 );
						else
							++it1;
					}
					path.front().second = 0;
				}
				r.push_back( std::move( path ) );
			}
			return r;
		}
		// topology_sketch_type
		template <
			typename CoordinateType, 
			typename PointType = boost::polygon::point_data<CoordinateType>
		>
		class topology_sketch_type 
		{
			// sub-types
			public:
				using coordinate_type					= CoordinateType;
				using point_type						= PointType;
				using linestring_type					= boost::geometry::model::linestring<point_type>;
			private:
				using voronoi_diagram_type				= boost::polygon::voronoi_diagram<CoordinateType>;
				using voronoi_cell_type					= typename voronoi_diagram_type::cell_type;
				using delaunay_triangulation_graph_type	= syc::topology::detail::v_01::delaunay_triangulation_graph_type<CoordinateType>;
				using sketchable_forest_type			= syc::topology::model::v_01::sketchable_forest<CoordinateType, PointType>;
				using net_type							= typename boost::graph_traits<typename sketchable_forest_type::net_type>::edge_descriptor; // very ungly!!
			// copy controls
			public:
				topology_sketch_type() = delete;	
				template <typename RootContainer, typename LeafContainer>
				topology_sketch_type( RootContainer&& roots, LeafContainer&& leafs ): 
					_num_roots( roots.size() ), _num_leafs( leafs.size() ), _SF( _num_roots + _num_leafs )
				{
					RootContainer	_roots( std::forward<RootContainer>( roots ) );
					LeafContainer	_leafs( std::forward<LeafContainer>( leafs ) );
					_points.insert( _points.end(), std::make_move_iterator( _roots.begin() ), std::make_move_iterator( _roots.end() ) );
					_points.insert( _points.end(), std::make_move_iterator( _leafs.begin() ), std::make_move_iterator( _leafs.end() ) );
					construct_delaunay_trianguation();
					construct_sketchable_forest();
				};
			// data members
			private:
				std::vector<point_type>					_points;
				std::size_t								_num_roots;
				std::size_t								_num_leafs;

				voronoi_diagram_type					_VD;
				delaunay_triangulation_graph_type		_DTG;
				sketchable_forest_type					_SF;			// For now, net is embeded in the forest. UGLY SYNTAX!!
			// method
			public:
				template <typename NetIter>
				unsigned the_basic_routing_algorithm( NetIter it, NetIter it_end );

				static 
				bool clockwise_depth_first_route( typename sketchable_forest_type::topo_descriptor& h, std::size_t const& u, sketchable_forest_type& F );

				template <typename NetIter>
				unsigned t_escape( NetIter it, NetIter it_end );

				bool search( typename sketchable_forest_type::topo_descriptor& h, std::size_t const& u, sketchable_forest_type& F, unsigned x );

				bool search2( bool direction, typename sketchable_forest_type::topo_descriptor& h, net_type e, unsigned x );

				std::vector<linestring_type> realization( void );

				template <typename string, typename VecLinestring> inline
				void plot( string name, VecLinestring const& GPs ) const; 

				// debug
				void print_capacities( void ) const
				{
					for ( auto ei = boost::edges( _DTG ); ei.first != ei.second; ++ei.first )
					{
						auto e = *ei.first;
						auto x = boost::get( boost::edge_capacity, _DTG, e );
						std::cout <<"("<<boost::source( e, _DTG )<<","<<boost::target( e, _DTG )<<",{"<<x<<"})"<< std::endl;
					}
				}
				void print_capacities( std::size_t u, std::size_t v ) const
				{
					auto e = boost::edge( u, v, _DTG );
					auto x = boost::get( boost::edge_capacity, _DTG, e );
					std::cout << "(" << u << "," << v << ",{" << x << "})" << std::endl;
				}
			private:
				inline 
				void construct_delaunay_trianguation( void );	// setup _VD, _TG

				inline 
				void construct_sketchable_forest( void );		// setup _SF based on _TG
		};
		// topology_sketch_type construct_delaunay_trianguation
		template <typename CoordinateType, typename PointType>
		inline void topology_sketch_type<CoordinateType, PointType>::
		construct_delaunay_trianguation( void )  
		{
			boost::polygon::construct_voronoi( _points.begin(), _points.end(), &_VD );
			for ( auto cell_it = _VD.cells().begin(); cell_it != _VD.cells().end(); ++cell_it ) {
				cell_it->color( 1 );
				auto	u		= cell_it->source_index();
				auto	edge_it	= cell_it->incident_edge();
				do {
					if ( edge_it->is_primary() ) {
						auto adj_cell_it	= edge_it->twin()->cell(); 
						if ( adj_cell_it->color() == 0 ) {
							auto v = adj_cell_it->source_index();
							auto e = boost::add_edge( u, v, _DTG ).first;
							if		( u < _num_roots && v < _num_roots ) {
								boost::put( boost::edge_weight, _DTG, e, std::numeric_limits<CoordinateType>::max() );
							}
							else {
								boost::put( boost::edge_weight, _DTG, e, boost::geometry::comparable_distance( _points[ u ], _points[ v ] ) );
							}
							boost::put( boost::edge_capacity, _DTG, e, 0 );
						}
					}
					edge_it = edge_it->prev();
				} while ( edge_it != cell_it->incident_edge() );
			}
			for ( auto cell_it = _VD.cells().begin(); cell_it != _VD.cells().end(); ++cell_it ) {
				boost::put( vertex_cell, _DTG, cell_it->source_index(), cell_it );
			}
		}
		// topology_sketch_type construct_sketchable_forest
		template <typename CoordinateType, typename PointType>
		inline void topology_sketch_type<CoordinateType, PointType>::
		construct_sketchable_forest( void ) 
		{
			// add virtual root
			std::size_t virtual_root = boost::add_vertex( _DTG );
			for ( std::size_t i = 0; i != _num_roots; ++i ) {
				boost::put( boost::edge_weight, _DTG, boost::add_edge( virtual_root, i, _DTG ).first, 0.0 );
			}
			// find edges of MST by kruskal algorithm and mark them in _TG
			std::vector<typename boost::graph_traits<delaunay_triangulation_graph_type>::edge_descriptor> 
				spanning_tree_edges;
			boost::kruskal_minimum_spanning_tree( _DTG, std::back_inserter( spanning_tree_edges ) );
			auto const EIDX = boost::get( boost::edge_index, _DTG );
			for ( auto const& e: spanning_tree_edges ) EIDX[ e ] = 1;		// index of 1 marks the edge as an edge of the MST
			// filter the MST out of _TG
			using EdgeIndxMap	= typename boost::property_map<delaunay_triangulation_graph_type, boost::edge_index_t>::type;
			using is_MST_edge	= typename syc::topology::detail::v_01::is_MST_edge<EdgeIndxMap>;
			boost::filtered_graph<delaunay_triangulation_graph_type, is_MST_edge> spannig_tree( _DTG, is_MST_edge( EIDX ) );
			// remove the virtual root in _TG 
			boost::clear_vertex( virtual_root, _DTG );
			boost::remove_vertex( virtual_root, _DTG );

			// adapt the undirected spanning tree to the sketchable-forest _SF which is oredered and directed
			for ( std::size_t i = 0; i != _num_roots; ++i ) {
				// breadth_first_search	the MST, a filtered graph of TS
				boost::breadth_first_search( spannig_tree, i , visitor( syc::topology::detail::v_01::build_sketchable_forest( _DTG, _SF ) ) );
			}
			_SF.make_sketch();
		}
		// topology_sketch_type the_basic_routing_algorithm
		template <typename CoordinateType, typename PointType>
		template <typename NetIter>
		unsigned topology_sketch_type<CoordinateType, PointType>::
		the_basic_routing_algorithm( NetIter it, NetIter it_end )
		{
			std::size_t failure_cnt = 0; 
			// add_two_pin_nets
			for ( ; it != it_end; ++it ) 
			{
				boost::add_edge( it->first, it->second, _SF.nets );
			}
			unsigned loop = 0;
			for ( auto ni = boost::edges( _SF.nets ); ni.first != ni.second; ++ni.first, ++loop )	// For all nets.
			{
				//if ( loop > 1 ) continue;
				auto source_vertex = boost::source( *ni.first, _SF.nets );					// Get the source vertex 
				auto target_vertex = boost::target( *ni.first, _SF.nets );					// and target vertex.	
				auto s = _SF.topo_instances( source_vertex )[ 0 ];							// Get the first topo_descriptor coorsponding to the source vertex.
				_SF.reset_slice_indexes();													// All the slices are not yet visited.
				if ( !clockwise_depth_first_route( s, target_vertex, _SF ) ) { ++failure_cnt; } 
			}
			return failure_cnt;
		}
		// topology_sketch_type clockwise_depth_first_route
		template <typename CoordinateType, typename PointType>
		bool topology_sketch_type<CoordinateType, PointType>::
		clockwise_depth_first_route( typename sketchable_forest_type::topo_descriptor& h, std::size_t const& u, sketchable_forest_type& F ) 
		{
			// Clockwise seach for the first target vertex t of h that is within the same slice of h.
			// If sucessfully finds a t, route from h to t and move on to the next source vertex.
			auto m = F.slice( h ); 																		// Let m = the pointer to slice where h is in.
			m->index = 1; 																				// Set the color of slice m gray
			for ( auto const& t : F.topo_instances( u ) ) { 
				auto n = F.slice( t ); 																	// Let n = the pointer to slice where target pointed by t is in.
				if ( m == n ) {
					F.make_slice( h, t );	
					return true; 																		// Successfully find a topological routing sketch
				}
			}
			// EVENT POINT: cannot route from h to any of its targets in the current slice
			
			// Anticlockwise search for any candidate edges to go through
			std::stack<typename sketchable_forest_type::topo_descriptor> E;
			auto ori_mark = F.slice( h )->end();
			for ( auto p = h; ++p != ori_mark; ) {
				if ( F.attribution( p ) != 2 ) { 
					auto e = F.ordering_pair( p );															// Let e be the related edge of the edge pointed by p
					if ( 0 == F.slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
						E.push( p );																		// the edge pointed by pi is chosen as a candidate.
					}
				}
			}
			for ( auto p = ori_mark; ++p != h; ) { 
				if ( F.attribution( p ) != 2 ) { 
					auto e = F.ordering_pair( p );															// Let e be the related edge of the edge pointed by p
					if ( 0 == F.slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
						E.push( p );																		// the edge pointed by pi is chosen as a candidate.
					}
				}
			}
			while ( !E.empty() ) {
				auto inArrow = E.top();
				h = F.make_slice( h, inArrow );					
				if ( clockwise_depth_first_route( h, u, F ) ) { 
					return true; 
				}
				E.pop();
			}

			// EVENT POINT: In the current slice, there is no target to route to, and no edge go through. I.e Have to turn around, back to the previous slice.
			
			h = F.free_slice( h );

			return false;
		}
		// topology_sketch_type t_escape
		template <typename CoordinateType, typename PointType>
		template <typename NetIter>
		unsigned topology_sketch_type<CoordinateType, PointType>::
		t_escape( NetIter it, NetIter it_end )
		{
			unsigned num_failed_nets = 0;
			// add_two_pin_nets
			for ( ; it != it_end; ++it ) 
			{
				boost::add_edge( it->first, it->second, _SF.nets );
			}
			unsigned loop = 0;
			for ( auto ni = boost::edges( _SF.nets ); ni.first != ni.second; ++ni.first, ++loop )	// For all nets.
			{	
				//if ( loop > 1 ) continue;
				auto e = *ni.first;
				auto source_vertex = boost::source( e, _SF.nets );							// Get the source vertex 
				auto target_vertex = boost::target( e, _SF.nets );							// Get the target vertex 
				auto x = source_vertex;
				while (boost::in_degree( x, _SF ) != 0) { 
					auto [in_edge_iterator_begin, in_edge_iterator_end] = boost::in_edges( x, _SF );
					x = boost::source( *in_edge_iterator_begin, _SF );
				} 
				auto head = _SF.topo_instances( source_vertex ).front();
				if ( !search( head, target_vertex, _SF, x ) ) 
				//if ( !search2( 0, head, e, x ) ) 
				{
					_SF.reset_slice_indexes();													// All the slices are not yet visited.
					if ( !clockwise_depth_first_route( head, target_vertex, _SF ) )
						++num_failed_nets;
				}
			}
			return num_failed_nets;
		}
		// topology_sketch_type search 
		template <typename CoordinateType, typename PointType>
		bool topology_sketch_type<CoordinateType, PointType>::
		search( typename sketchable_forest_type::topo_descriptor& h, std::size_t const& u, sketchable_forest_type& F, unsigned x )
		{
			F.txt_vitualizer();
			auto f = [&]( unsigned r ) -> unsigned { return (_num_roots / 2 - x + r + _num_roots) % _num_roots; };
			auto m 		= F.slice( h );
			std::stack<typename sketchable_forest_type::topo_descriptor> E;
			if ( 1 == F.attribution( h ) ) {		// a back side of an edge
				// searching dounward 
				for ( auto scout = h; ++scout != m->end(); ) { 
					auto attr = F.attribution( scout );
					if 		( 0 == attr ) { 		// It's about to passed by the back side of an edge. I.e. it is about to search in a leaf-oriented way.
						h = F.make_slice( h, scout );
						if ( search( h, u, F, x ) ) return true;
					}
					else if ( 1 == attr ) {
						E.push( scout );	
					}
					else if ( 2 == attr ) {
						auto v = F.corresponding_vertex( scout );
						if ( 0 != boost::in_degree( v, F ) ) {	// it is not a root
							v = *boost::adjacent_vertices( v, F.nets ).first;
						}
						auto fv = f(v);
						auto fu = f(u);
						if ( f(v) > f(u) ) {
							if ( std::prev( scout ) != h ) {
								h = F.make_slice( h, std::prev( scout ) );
								if ( search( h, u, F, x ) ) return true;
							}
							h = F.free_slice( h );
							return false;
						}
						else if ( v == u )	{
							h = F.make_slice( h, scout );
							return true;
						}
					}
				}
				for ( auto scout = h; !E.empty(); E.pop() ) {
					scout = E.top();
					h = F.make_slice( h, scout );
					if ( search( h, u, F, x ) ) return true;
				}
			}
			else {									// a vertex or a tree side of an edge
				// searching upward 
				for ( auto scout = h; --scout != m->end(); ) { 
					auto attr = F.attribution( scout );
					if 		( 0 == attr ) {
						E.push( scout );	
					}
					else if ( 1 == attr ) { 		// It's about to passed by the back side of an edge. I.e. it is about to search in a leaf-oriented way.
						h = F.make_slice( h, scout );
						if ( search( h, u, F, x ) ) return true;
					}
					else if ( 2 == attr ) {
						auto v = F.corresponding_vertex( scout );
						if ( 0 != boost::in_degree( v, F ) ) {	// it is not a root
							v = *boost::adjacent_vertices( v, F.nets ).first;
						}
						auto fv = f(v);
						auto fu = f(u);
						if ( f(v) < f(u) ) {
							if ( std::next( scout ) != h ) {
								h = F.make_slice( h, std::next( scout ) );
								if ( search( h, u, F, x ) ) return true;
							}
							h = F.free_slice( h );
							return false;
						}
						else if ( v == u )	{
							h = F.make_slice( h, scout );
							return true;
						}
					}
				}
				for ( auto scout = h; !E.empty(); E.pop() ) {
					scout = E.top();
					h = F.make_slice( h, scout );
					if ( search( h, u, F, x ) ) return true;
				}
			}
			h = F.free_slice( h );
			return false;
		}
		// topology_sketch_type search2 
		template <
			typename CoordinateType, 
			typename PointType
		>
		bool 
		topology_sketch_type<CoordinateType, PointType>::search2( 
			bool direction, typename sketchable_forest_type::topo_descriptor& h, net_type e, unsigned x 
		) {

			auto f = [&]( unsigned r ) -> unsigned { return (_num_roots / 2 - x + r + _num_roots) % _num_roots; };
			_SF.txt_vitualizer();
			auto u = boost::target( e, _SF.nets );
			// searching dounward (counterclockwise) 
			if (auto m = _SF.slice( h ); direction == 0) {								 	
				for ( auto scout = h; ++scout != h; ) { 
					if (scout == m->end()) { continue; }

					auto attr = _SF.attribution( scout );
					if 		( 0 == attr ) {						// It's about to passed by the back side of an edge. I.e. it is about to search in a leaf-oriented way.
						h = _SF.make_slice( h, scout );
						if ( search2( 0, h, e, x ) ) return true;
					}
					else if ( 1 == attr ) {
						;
					}
					else if ( 2 == attr ) {
						auto v = _SF.corresponding_vertex( scout );
						if ( v == u )	{
							h = _SF.make_slice( h, scout );
							return true;
						}
						else if ( 0 != boost::in_degree( v, _SF ) ) {				// it is not a root
							v = *boost::adjacent_vertices( v, _SF.nets ).first;
						}

						if (f( v ) > f( u )) {
							if ( std::prev( scout ) != h ) {
								h = _SF.make_slice( h, std::prev( scout ) );
								if ( search2( 1, h, e, x ) ) return true;
							}
							h = _SF.free_slice( h );
							return false;
						}
					}
				} 
			}
			// searching upward (clockwise) 
			else {									
				for ( auto scout = h; --scout != m->end(); ) { 
					if (scout == m->end()) { continue; }

					auto attr = _SF.attribution( scout );
					if 		( 0 == attr ) {
						;
					}
					else if ( 1 == attr ) { 		// It's about to passed by the back side of an edge. I.e. it is about to search in a leaf-oriented way.
						h = _SF.make_slice( h, scout );
						if ( search2( 1, h, e, x ) ) return true;
					}
					else if ( 2 == attr ) {
						auto v = _SF.corresponding_vertex( scout );
						if ( v == u )	{
							h = _SF.make_slice( h, scout );
							return true;
						}
						else if ( 0 != boost::in_degree( v, _SF ) ) {	// it is not a root
							v = *boost::adjacent_vertices( v, _SF.nets ).first;
						}

						if (f( v ) < f( u )) {
							if ( std::next( scout ) != h ) {
								h = _SF.make_slice( h, std::next( scout ) );
								if ( search2( 0, h, e, x ) ) return true;
							}
							h = _SF.free_slice( h );
							return false;
						}
					}
				}
			}
			h = _SF.free_slice( h );
			return false;
		}
		// topology_sketch_type realization 
		template <typename CoordinateType, typename PointType>
		inline auto topology_sketch_type<CoordinateType, PointType>::
		realization( void ) -> std::vector<linestring_type>
		{
			auto R = _SF.results();															// every net should have a head before calling this method
			std::vector<std::vector<std::tuple<unsigned, unsigned, int>>>	TCPs;			// triangulation crossing paths
			// for all paths in the results R from _SF
			for ( auto& path : R ) 
			{						
				// step1: relax each of the paths
				for ( bool relax = true; relax == true; ) {
					auto it1 = path.begin();
					relax = false;
					while ( true ) {
						auto it2 = std::next( it1 );
						if ( it2 == path.end() )	break;
						if ( it2->second == 0 )		break;
						if ( it2->second > 0 ) {
							auto cell_it = boost::get( vertex_cell, _DTG, it1->first );
							auto vor_edge_it = cell_it->incident_edge();
							while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->next();	
							for ( ;    std::next( it2 )->first == vor_edge_it->next()->twin()->cell()->source_index()
									&& std::next( it2 )->second >= 0;
								  vor_edge_it = vor_edge_it->next() ) {
								it2 = path.erase( it2 );
								relax = true;
							}
						}
						else if ( it2->second < 0 ) {
							auto cell_it = boost::get( vertex_cell, _DTG, it1->first );
							auto vor_edge_it = cell_it->incident_edge();
							while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->prev();
							for ( ;    std::next( it2 )->first == vor_edge_it->prev()->twin()->cell()->source_index() 
									&& std::next( it2 )->second <= 0;
								  vor_edge_it = vor_edge_it->prev() ) {
								it2 = path.erase( it2 );
								relax = true;
							}
						}
						it1 = it2;
					}
				}
				// step2: realize each of the paths to a triangulation crossing paths
				std::vector<std::tuple<unsigned, unsigned, int>> X;
				if ( path.empty() )
				{
					TCPs.push_back( std::move( X ) );
				} 
				else
				{
					X.emplace_back( path.front().first, path.front().first, 0 );
					auto	it1			= path.begin(), 
							it2			= std::next( it1 ),
							it3			= std::next( it2 );
					auto	vor_edge_it	= boost::get( vertex_cell, _DTG, it1->first )->incident_edge();
					while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) 
					{
						vor_edge_it = vor_edge_it->next();
					}
					auto cap = boost::get( boost::edge_capacity, _DTG );
					for	( ;it2->second != 0; ++it1, ++it2, ++it3 )
					{
						vor_edge_it = vor_edge_it->twin();
						if ( it2->second > 0 ) 
						{ 
							if ( it1->second < 0 )
							{
								X.emplace_back( it2->first, it1->first, it2->second );
								auto x = ++cap[ boost::edge( it2->first, it1->first, _DTG ).first ];
							}
							vor_edge_it = vor_edge_it->prev();
							for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
								  v != it3->first; 
								  vor_edge_it = vor_edge_it->prev(), v = vor_edge_it->twin()->cell()->source_index() )  
							{
								X.emplace_back( it2->first, v, it2->second );
								auto x = ++cap[ boost::edge( it2->first, v, _DTG ).first ];
							}
						} 
						else if ( it2->second < 0 ) 
						{
							if ( it1->second > 0 )
							{
								X.emplace_back( it2->first, it1->first, it2->second );
								auto x = ++cap[ boost::edge( it2->first, it1->first, _DTG ).first ];
							}
							vor_edge_it = vor_edge_it->next();
							for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
								  v != it3->first; 
								  vor_edge_it = vor_edge_it->next(), v = vor_edge_it->twin()->cell()->source_index() )  
							{
								X.emplace_back( it2->first, v, it2->second );
								auto x = ++cap[ boost::edge( it2->first, v, _DTG ).first ];
							}
						}
					}
					X.emplace_back( path.back().first, path.back().first, 0 );
					TCPs.push_back( std::move( X ) );
				}
			}
			unsigned loop = 0;
			std::vector<linestring_type> GPs;			// geometrical paths
			for ( auto const& X : TCPs ) 
			{		
				// step3: realize each of the triangulation crossing paths into a geometrical path
				linestring_type gpath;
				for ( auto const& t : X ) { 
					auto const& u = std::get<0>(t);
					auto const& v = std::get<1>(t);
					auto const& w = std::get<2>(t);
					auto const& p0 = _points[ u ];
					auto const& p1 = _points[ v ];
					auto ei = boost::edge( u, v, _DTG );
					if ( u == v )
					{
						auto const& x = boost::geometry::get<0>(p0);
						auto const& y = boost::geometry::get<1>(p0);
						boost::geometry::append( gpath, point_type( x, y ) );
					}
					else 
					{
						double b = boost::get( boost::edge_capacity, _DTG, boost::edge( u, v, _DTG ).first ) + 1.5;			// 1 beacause of intervals, 0.5 for reservation for special cases
						auto frac = static_cast<double>(std::abs(w)) / b; 
						auto const& x = boost::geometry::get<0>(p0) - ( boost::geometry::get<0>(p0) - boost::geometry::get<0>(p1) ) * frac;
						auto const& y = boost::geometry::get<1>(p0) - ( boost::geometry::get<1>(p0) - boost::geometry::get<1>(p1) ) * frac;
						boost::geometry::append( gpath, point_type( x, y ) );
					}
				}
				GPs.push_back( std::move( gpath ) );
				++loop;
			}
			/*
			// ==============================debug_start========================================
			std::cout << "R" << std::endl;
			std::cout << std::string( 20, '-' ) << std::endl;
			for ( unsigned i = 0; i != R.size(); ++i )
			{
				if ( i == 0 ) 
				{ 
					std::cout << "[" << i << "]" << std::endl;
					for ( auto const& x: R[i] )
					{
						std::cout << "(" << x.first << "," << x.second << ") ";	
					}
					std::cout << std::endl;
				}
			}
			std::cout << "relaxed R" << std::endl;
			std::cout << std::string( 20, '-' ) << std::endl;
			for ( unsigned i = 0; i != R.size(); ++i )
			{
				if ( i == 0 ) 
				{ 
					std::cout << "[" << i << "]" << std::endl;
					for ( auto const& x: R[i] )
					{
						std::cout << "(" << x.first << "," << x.second << ") ";	
					}
					std::cout << std::endl;
				}
			}
			std::cout << "TCPs" << std::endl;
			std::cout << std::string( 20, '-' ) << std::endl;
			for ( unsigned i = 0; i != R.size(); ++i )
			{
				if ( i == 0 ) 
				{ 
					std::cout << "[" << i << "]" << std::endl;
					for ( auto const& x: TCPs[i] )
					{
						std::cout << "(" << std::get<0>( x ) << "," << std::get<1>( x ) << "," << std::get<2>( x ) << ") ";	
					}
					std::cout << std::endl;
				}
			}
			// ==============================debug_end========================================
			*/
			return GPs;
		}
		// topology_sketch_type plot 
		template <typename CoordinateType, typename PointType>
		template <typename string, typename VecLinestring>
		inline void topology_sketch_type<CoordinateType, PointType>::
		plot( string name, VecLinestring const& GPs ) const 
		{
			// vitualization using gnugraph
			syc::gnuplot::v_01::gnugraph G( name );
			G.set_title( name );
			// create data files
			std::shared_ptr<std::ofstream> out;
			// write data
			out = G.plot_datafile_with_lines( "./tri_segs.dat", "", G.set_style_line( "lt 7 lc rgb 'red' ps 10" ) );
			for ( auto range = boost::edges( _DTG ); range.first != range.second; ++range.first ) {
				auto e = *range.first;
				auto const& source_point = _points[ boost::source( e, _DTG ) ];
				auto const& target_point = _points[ boost::target( e, _DTG ) ];
				*out << boost::geometry::get<0>( source_point ) << " " << boost::geometry::get<1>( source_point ) << "\n"
					 << boost::geometry::get<0>( target_point ) << " " << boost::geometry::get<1>( target_point ) << "\n\n";
			}
			out->close();
			out = G.plot_datafile_with_lines( "./tree_segs.dat", "", G.set_style_line( "lt 7 lc rgb 'blue' ps 10" ) );
			for ( auto range = boost::edges( _SF ); range.first != range.second; ++range.first ) {
				auto e = *range.first;
				auto const& source_point = _points[ boost::source( e, _SF ) ];
				auto const& target_point = _points[ boost::target( e, _SF ) ];
				*out << boost::geometry::get<0>( source_point ) << " " << boost::geometry::get<1>( source_point ) << "\n"
					 << boost::geometry::get<0>( target_point ) << " " << boost::geometry::get<1>( target_point ) << "\n\n";
			}
			out->close();
			out = G.plot_datafile_with_labels_point( "./points.dat", "", G.set_style_line( "lt 7 lc rgb 'gray'" ) );
			for ( auto range = boost::vertices( _SF ); range.first != range.second; ++range.first ) {
				auto v = *range.first;
				auto const& p = _points[ v ];
				*out << boost::geometry::get<0>( p ) << " " << boost::geometry::get<1>( p ) << " " << v << "\n";
			}
			out->close();
			unsigned loop = 0;
			out = G.plot_datafile_with_lines( "./paths.dat", "", G.set_style_line( "lt 7 lc rgb 'green' ps 10" ) );
			for ( auto const& gp : GPs ) 
			{		 
				for ( auto const& point : gp ) 
				{
					*out << boost::geometry::get<0>(point) << " " << boost::geometry::get<1>(point) << "\n";
				}
				*out << "\n";
			}
			out->close();
			// close scripts and system call to gnuplot
			G.plot( "-p" );
			// clean all the scripts and data files
			G.clean_all();
		}
	}
}
namespace syc::topology::v_02 
{
	using namespace fmt::literals;

	// # design detail
	
	namespace detail 
	{
		// ## predicate: is_triangulation_edge
		
		template <typename EdgeIndexMap>
		struct is_triangulation_edge 
		{
			public:
				is_triangulation_edge() = default;
				is_triangulation_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 0;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ## predicate: is_MST_edge
		
		template <typename EdgeIndexMap>
		struct is_MST_edge 
		{
			public:
				is_MST_edge() = default;
				is_MST_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 1;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ## predicate: is_boundary_edge
		
		template <typename EdgeIndexMap>
		struct is_boundary_edge 
		{
			public:
				is_boundary_edge() = default;
				is_boundary_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 2;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ## predicate: is_net_edge
		
		template <typename EdgeIndexMap>
		struct is_net_edge 
		{
			public:
				is_net_edge() = default;
				is_net_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 3;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ## set_parent_vertex_property in property map of vertex_mother_t
		
		template <typename PMap>
		struct set_mother_vertex_property: public boost::default_bfs_visitor 
		{
			// copy control
			public:
				set_mother_vertex_property( PMap m ): _m( m ) {}
			// data member
			private:
				PMap _m;
			// events
			public:
				// tree_edge
				template <typename Edge, typename Graph>	
				void tree_edge( Edge e, Graph const& g )  
				{
					auto s = source( e, g );
					auto t = target( e, g );
					_m[ t ] = s;
				}
		};
		
		// ## accumulate_cross_number
		
		template < typename Net >
		struct accumulate_cross_number : public boost::default_dfs_visitor
		{
			// copy controls
				accumulate_cross_number( std::size_t const r, std::size_t const m, Net const& N, unsigned& n ): 
					_num_roots( r ),
					_nets( N ),
					_seq( m ),
					_cross_number( n )
				{}

			// data members
				std::size_t									_num_roots;			// IN
				Net const&									_nets;				// IN
				std::vector< std::list< std::size_t > >		_seq;				// UTIL 
				std::size_t									_current_root;		// UTIL
				unsigned&									_cross_number;		// OUT 

            // start_vertex
            template < typename Vertex, typename Graph >
            void start_vertex( Vertex u, Graph const& g ) 
            {
                _current_root = u;
            }

            // finish_vertex
            template < typename Vertex, typename Graph >
            void finish_vertex( Vertex u, Graph const& g ) 
            {
                auto dis = [&]( auto v ) -> int
                {
                    int Dis_CW = ( _num_roots + _current_root - *adjacent_vertices( v, _nets ).first ) % _num_roots;
                    int Dis_CCW = _num_roots - Dis_CW;
                    return ( Dis_CW < Dis_CCW ) ? Dis_CW : - Dis_CCW;
                };
                auto comp = [&]( auto f2, auto f1 ) 
                {
                    if ( dis(f2) < dis(f1) ) 
                    {
                        ++_cross_number;
                        return true;
                    }
                    else 
                    {
                        return false;
                    }
                };
                auto comp2 = [&]( auto f2, auto f1 ) 
                {
                    return ( dis(f2) < dis(f1) );
                };
                auto merge_sort_and_compute_inversion = [&]( auto w, auto rep )
                {
                    if ( rep == 0 ) return;

                    auto edge_it = get( vertex_cell, g, u )->incident_edge();
                    for ( 
                        auto v = edge_it->twin()->cell()->source_index(); 
                        v != w; 
                        edge_it = edge_it->prev(), v = edge_it->twin()->cell()->source_index()
                    ) ;

                    for ( int i = 0; i != rep; ) 
                    {
                        edge_it = edge_it->prev();
                        auto v = edge_it->twin()->cell()->source_index();
                        if ( edge( u, v, g ).second ) 
                        {
                            auto middle = _seq[ u ].begin();
                            _seq[u].splice( _seq[u].begin(), _seq[v] );
                            std::inplace_merge( _seq[u].begin(), middle, _seq[u].end(), comp ); 			
                            ++i;
                        }
                    } 
                };
                auto push_and_sort = [&] ()
                {
                    auto middle = _seq[ u ].begin();
                    _seq[u].push_front( u );
                    std::inplace_merge( _seq[u].begin(), middle, _seq[u].end(), comp2 ); 			
                };
                auto w = get( vertex_mother, g, u );
                auto deg = out_degree( u, g );
                if ( w == u )				// it's a root 
                {
                    merge_sort_and_compute_inversion( (w - 1) % _num_roots, deg );

                    auto l = dis( _seq[u].front() ), r = dis( _seq[u].back() ); 
                    auto sz = _seq[u].size();
                    auto m = r - l + 1 - sz;
                    _cross_number += m * (sz + 1) / 2;
                }
                else if ( deg == 1 )					// it's a leaf
                {
                    _seq[ u ].emplace_back( u );
                }
                else												// it's the others
                {
                    merge_sort_and_compute_inversion( w, deg - 1 );
                    push_and_sort();
                }
            }
		};
	
        // ## compute_dest_ordering_and_cross_number

		template < typename Net >
        struct compute_dest_ordering_and_cross_number : public boost::default_dfs_visitor 
        {
            // declarition
            using DestOrder = std::vector< std::pair< int, std::size_t > >;
            using SortingStack = std::vector< std::list< std::size_t > >;
            // data member
            Net const&              N;          // IN:      nets
            int                     B;          // IN:      number of roots on the boudary
            int                     R;          // UTIL:    current root
            DestOrder               D;          // UTIL:    dest_order
            SortingStack            S;          // UTIL:    sort to compute cross number
            unsigned&               cross_n;    // OUT:     cross_number
            // constructor
            compute_dest_ordering_and_cross_number( Net const& n, int num_vertices, int num_roots, unsigned& cn ):
                N( n ), B( num_roots ), D( num_vertices ), cross_n( cn ) {}
            // start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
            // finish_vertex
			template < typename Vertex, typename Graph >
			void finish_vertex( Vertex u, Graph const& g ) 
			{
                compute_D( u );

                if ( out_degree( u, g ) == 0  ) {
                    if ( u >= B ) { S.emplace_back( 1, u ); }   // it is a leaf
                }
                else
                {
                    auto comp = [&]( auto f2, auto f1 ) 
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        return d2 < d1;
                    };
                    unsigned f1_left;
                    auto comp2 = [&]( auto f2, auto f1 ) mutable
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        if ( d2 < d1 )
                        {
                            cross_n += f1_left;
                            return true;
                        }
                        else
                        {
                            --f1_left;
                            return false;
                        }
                    };
                    std::list<std::size_t> L( std::move( S.back() ) );
                    S.pop_back(); 
                    auto merge_and_compute_cross_n = [&] ( void )
                    {
                        f1_left = L.size();
                        L.merge( std::move( S.back() ), comp2 );
                        S.pop_back(); 
                    };
                    for ( unsigned i = 0, i_end = out_degree( u, g ) - 1; i != i_end; ++i ) 
                    {
                        merge_and_compute_cross_n();
                    }
                    if ( u < B ) 
                    {
                        if ( !S.empty() ) 
                        {
                            merge_and_compute_cross_n();
                        }
                    }
                    else 
                    {
                        L.merge( std::list<std::size_t>( 1, u ), comp );
                    }
                    S.push_back( std::move( L ) );
                }
                if ( u == B - 1 )
                {
                    unsigned f2_size;
                    auto comp3 = [&]( auto f2, auto f1 ) mutable
                    {
                        if ( f2 > f1 )
                        {
                            --f2_size;
                            return true;
                        }
                        else
                        {
                            cross_n += f2_size;
                            return false;
                        }
                    };
                    std::list<std::size_t> X;
                    for ( auto const& v: S.back() ) { X.emplace_back( *adjacent_vertices( v, N ).first ); }
                    auto offset = B - X.back();
                    for ( auto& v: X ) { v = (v + offset) % B; }
                    auto it = std::next( std::find( X.begin(), X.end(), 1 ) );
                    f2_size = X.size() - std::distance( X.begin(), it );
                    std::inplace_merge( X.begin(), it, X.end(), comp3 );
                }
            }
            // compute_D
			template < typename Vertex >
            void compute_D( Vertex u )
            {
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    std::get<0>( D[ u ] ) = ( Dcw < Dccw )? Dcw : -Dccw;
                    std::get<1>( D[ u ] ) = R;
                }
                else
                {
                    std::get<0>( D[ u ] ) = 0;
                    std::get<1>( D[ u ] ) = R;
                } 
            }
        };

        // ## compute_dest_ordering_and_cross_number2

		template < typename Net >
        struct compute_dest_ordering_and_cross_number2 : public boost::default_dfs_visitor 
        {
            // declarition
            using DestOrder = std::vector< std::pair< int, std::size_t > >;
            using SortingStack = std::vector< std::list< std::size_t > >;
            // data member
            Net const&              N;          // IN:      nets
            int                     B;          // IN:      number of roots on the boudary
            int                     R;          // UTIL:    current root
            DestOrder               D;          // UTIL:    dest_order
            SortingStack            S;          // UTIL:    sort to compute cross number
            unsigned&               cross_n;    // OUT:     cross_number
            // constructor
            compute_dest_ordering_and_cross_number2( Net const& n, int num_vertices, int num_roots, unsigned& cn ):
                N( n ), B( num_roots ), D( num_vertices ), cross_n( cn ) {}
            // start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
            // finish_vertex
			template < typename Vertex, typename Graph >
			void finish_vertex( Vertex u, Graph const& g ) 
			{
                compute_D( u );

                if ( out_degree( u, g ) == 0  ) {
                    if ( u >= B ) { S.emplace_back( 1, u ); }   // it is a leaf
                }
                else
                {
                    auto comp = [&]( auto f2, auto f1 ) 
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        return d2 < d1;
                    };
                    unsigned f1_left;
                    auto comp2 = [&]( auto f2, auto f1 )
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        if ( d2 < d1 )
                        {
                            cross_n += f1_left;
                            return true;
                        }
                        else
                        {
                            --f1_left;
                            return false;
                        }
                    };
                    std::list<std::size_t> L( std::move( S.back() ) );
                    S.pop_back(); 
                    auto merge_and_compute_cross_n = [&] ( void )
                    {
                        f1_left = L.size();
                        L.merge( std::move( S.back() ), comp2 );
                        S.pop_back(); 
                    };
                    for ( unsigned i = 0, i_end = out_degree( u, g ) - 1; i != i_end; ++i ) 
                    {
                        merge_and_compute_cross_n();
                    }
                    if ( u < B ) 
                    {
                        auto l = std::get<0>( D[ L.front() ] ) - long( std::get<1>( D[ L.front() ] ) );
                        auto r = std::get<0>( D[ L.back() ] ) - long( std::get<1>( D[ L.back() ] ) );
                        auto sz = L.size();
                        auto m = r - l + 1 - sz;
                        cross_n += m * (sz + 1) / 2;
                    }
                    else 
                    {
                        L.merge( std::list<std::size_t>( 1, u ), comp );
                    }
                    S.push_back( std::move( L ) );
                }
            }
            // compute_D
			template < typename Vertex >
            void compute_D( Vertex u )
            {
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    std::get<0>( D[ u ] ) = ( Dcw < Dccw )? Dcw : -Dccw;
                    std::get<1>( D[ u ] ) = R;
                }
                else
                {
                    std::get<0>( D[ u ] ) = 0;
                    std::get<1>( D[ u ] ) = R;
                } 
            }
        };

		// ## build_sketchable_forest
		
		template < typename sketable_forest >
		class build_sketchable_forest: public boost::default_bfs_visitor 
		{
			// constructor 
			public:
				build_sketchable_forest( sketable_forest& sf ): _SF( sf ) {}
			// data member
			private:
				sketable_forest&    _SF;
			// events
			public:	
				// examine_vertex
				template <typename Vertex, typename Graph>	
					void examine_vertex( Vertex u, Graph const& g )  
					{
						// reset color of all adjacent edges
						for ( auto ei = boost::out_edges( u, g ); ei.first !=  ei.second; ++ei.first ) {
							boost::put( boost::edge_color, g, *ei.first, 0 );
						}
					}
				// tree_edge
				template <typename Edge, typename Graph>	
					void tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 1 );
					}
				// non_tree_edge
				template <typename Edge, typename Graph>	
					void non_tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 2 );
					}
				// finish_vertex
				template <typename Vertex, typename Graph>	
					void finish_vertex( Vertex u, Graph const& g )  
					{
						auto it_end	= boost::get( vertex_cell, g, u )->incident_edge();
						auto it		= it_end->prev();
						while ( it != it_end ) {
							auto v = it->twin()->cell()->source_index();
							auto e = boost::edge( u, v, g ).first;
							if ( 2 == boost::get( boost::edge_color, g, e ) ) break;
							it = it->prev();
						}
						auto it2 = it;
						do {
							auto v = it2->twin()->cell()->source_index();
							auto e = boost::edge( u, v, g ).first;
							if ( 1 == boost::get( boost::edge_color, g, e ) ) {
								boost::add_edge( u, v, _SF );
							}
							it2 = it2->prev();
						} while ( it2 != it );
					}
		};

        // ## determine_net_priority_and_copmute_dest_ordering
        
		template < typename CircularFrame, typename Net, typename DestOrder >
		struct determine_net_priority_and_copmute_dest_ordering: public boost::default_dfs_visitor 
		{
            CircularFrame&  C;
            Net             N;
            int             B;      // number of roots on the boudary
            int             R;      // current root
            DestOrder&      D;

            determine_net_priority_and_copmute_dest_ordering( CircularFrame& c, Net& n, unsigned b, DestOrder& d ):
                C( c ), 
                N( n ), 
                B( b ),
                D( d )
            {}

			// start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
			// discover_vertex
			template < typename Vertex, typename Graph >
			void discover_vertex( Vertex u, Graph const& g ) 
			{
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    add_edge( u, *beg, C.nets ); 
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    D[ u ][ 0 ] = ( Dcw < Dccw )? Dcw : -Dccw;
                    D[ u ][ 1 ] = R;
                }
                else
                {
                    D[ u ][ 0 ] = 0;
                    D[ u ][ 1 ] = R;
                } 
            }
        };
        
	}

	// # planar_orderd_forest_solution
	
	template <
		typename CoordinateType, 
		typename PointType			= boost::polygon::point_data<CoordinateType>,
		typename VoronoiDiagramType	= boost::polygon::voronoi_diagram<CoordinateType>	
	>
	class planar_orderd_forest_solution:
		// ## public inherit a boost undirected graph. It's associated properties is collected here.
		public 
			boost::adjacency_list<
				boost::vecS, 

				boost::vecS, 

				boost::undirectedS,

				boost::property<vertex_idx_t,			std::size_t,		                // used to indicate the parent of the vertex 
				boost::property<vertex_mother_t,		std::size_t,		                // used to indicate the parent of the vertex 
				boost::property<vertex_deg_t,			int,		                        
				boost::property<vertex_cell_t,			typename boost::polygon::voronoi_diagram<CoordinateType>::const_cell_iterator
				>>>>,

				boost::property<boost::edge_index_t,	int,				                // Used as a marker. 0: Triangulation Graph Edge, 1: MST Edge, 2: Boundary, 3: net
				boost::property<boost::edge_color_t,	int,				                // Used as a tempary marker for traversing
				boost::property<boost::edge_weight_t,	CoordinateType,		                // Used by kruskal's algorithm
				boost::property<boost::edge_flow_t,		int,			
				boost::property<boost::edge_capacity_t,	int,				                // Used as a recording of # of traces passing an edge
				boost::property<edge_point_t,	        std::array<std::list<int>, 2>		// Used as a recording of geometry passing points. The front side of the list is the lower id vertex
				>>>>>>
			> 
	{
		// ## sub-types
		public:
			
			using coordinate_type		= CoordinateType;
			using point_type			= PointType;
			using linestring_type		= boost::geometry::model::linestring<point_type>;
			using voronoi_diagram_type	= VoronoiDiagramType;
	
		// ## data members
		public:
			std::size_t			num_roots;
			std::size_t			num_leafs;
		private:
			std::vector<PointType>					_points;
			std::shared_ptr<VoronoiDiagramType>		_vd;

            std::vector< unsigned >                 _perturbs; 
			
		// ## function members 

        // ### constructor

		public:
			template <typename RootContainer, typename LeafContainer, typename NetContainer>
				planar_orderd_forest_solution( RootContainer&& roots, LeafContainer&& leafs, NetContainer&& nets );
		private:
			void construct_delaunay_trianguation( void );
			void mark_Planar_Ordered_Forest_edges( void );
			template < typename NetContainer >
				void mark_net_edges( NetContainer const& nets );

        // ### gnuplot

        public:
			template < typename GeometricPath = std::vector< linestring_type > >
				void gnuplot( GeometricPath const& gp = GeometricPath{} );
        private:

        // ### SA functions

        public:
			auto point( std::size_t u ) { return this->_points[ u ]; }
			void print( std::ostream& os = std::cout );
			void perturb( void );
			bool switch_mother( std::size_t u );
			double cost( void );

        // ### write_perturb
        
        public:
            void write_perturb( std::ostream& os ) { for ( auto const& i : _perturbs ) { os << i << "\n"; } }

	};                                        

	// ## Public Methods of planar_orderd_forest_solution
	
	// ### constructor 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template <typename RootContainer, typename LeafContainer, typename NetContainer>
	planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	planar_orderd_forest_solution( RootContainer&& roots, LeafContainer&& leafs, NetContainer&& nets ): 
		num_roots( roots.size() ), num_leafs( leafs.size() ) 
	{
		RootContainer	_roots( std::forward<RootContainer>( roots ) );
		LeafContainer	_leafs( std::forward<LeafContainer>( leafs ) );

		_points.insert( _points.end(), std::make_move_iterator( _roots.begin() ), std::make_move_iterator( _roots.end() ) );
		_points.insert( _points.end(), std::make_move_iterator( _leafs.begin() ), std::make_move_iterator( _leafs.end() ) );

		_vd = std::make_shared<VoronoiDiagramType>(); 

		boost::polygon::construct_voronoi( _points.begin(), _points.end(), &*_vd );

		construct_delaunay_trianguation();

		mark_Planar_Ordered_Forest_edges();
		
		mark_net_edges( nets );
	};

	// ### print
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	print( std::ostream& os )
	{
		auto const& g = *this;
		unsigned col_width[5] = { 7, 35, 15, 13, 28 };

		auto header = 
			fmt::format( 
				"|-{0:-<{6}}-|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|\n" 
				"| {1: <{6}} | {2: ^{7}} | {3: <{8}} | {4: <{9}} | {5: <{10}} |\n"
				"|-{0:-<{6}}-|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|\n", 

				"", "vertex", "adjacency_list", "vertex_mother_t", "vertex_deg_t", "point get from vertex_cell_t",  
				col_width[0], col_width[1],	col_width[2], col_width[3], col_width[4]
			);

		os << header;

		for ( auto [it, it_end] = vertices( g ); it != it_end; ++it ) {

			auto [adj_it, adj_it_end] = adjacent_vertices( *it, g );
			std::vector<unsigned> adj_v( adj_it, adj_it_end );
			auto a1 = fmt::format( "{0}", adj_v );
			
			auto const& p = g._points[ get( vertex_cell, g, *it )->source_index() ];
			auto a2 = fmt::format( "({0: >10.2f}, {1: >10.2f})", p.x(), p.y() );

			auto row = 
				fmt::format(
					"| {0: <{5}} | {1: ^{6}} | {2: >{7}} | {3: >{8}} | {4: >{9}} |", 

					*it,			a1,				get( vertex_mother, g, *it ),	get( vertex_deg, g, *it ),	a2,  
					col_width[0],	col_width[1],	col_width[2],					col_width[3],				col_width[4]
				);

			os << row;	

			if (*it == (g.num_roots)) {
				os << fmt::format( "  <- num_roots: {}, num_leafs: {}", g.num_roots, g.num_leafs );
			}

			os << "\n";

		}


		unsigned col_width2[6] = { 10, 12, 12, 13, 11, 15 };

		auto header2 = 
			fmt::format( 
				"|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|-{0:-<{11}}-|-{0:-<{12}}-|\n" 
				"| {1: >{7}} | {2: >{8}} | {3: >{9}} | {4: >{10}} | {5: >{11}} | {6: >{12}} |\n" 
				"|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|-{0:-<{11}}-|-{0:-<{12}}-|\n" 

				, "", "edge", "edge_index_t", "edge_color_t", "edge_weight_t", "edge_flow_t", "edge_capacity_t"  
				, col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
			);

		os << header2;

		for (auto [it, it_end] = edges( g ); it != it_end; ++it) {

			auto edge = 
				fmt::format( "{{{0: >3}, {1: >3}}}", source( *it, g ), target( *it, g ) );

			auto w = get( boost::edge_weight, g, *it );
			auto weight = (w == std::numeric_limits<CoordinateType>::max())?  
				fmt::format( "{0: >{1}}", "max", col_width2[3] ): fmt::format( "{0: >{1}.2f}", w, col_width2[3] );
					

			auto row = 
				fmt::format( 
					"| {1: >{7}} | {2: >{8}} | {3: >{9}} | {4} | {5: >{11}} | {6: >{12}} |\n"

					, "", edge, get( boost::edge_index, g, *it ), get( boost::edge_color, g, *it ), weight, get( boost::edge_flow, g, *it ), get( boost::edge_capacity, g, *it ) 
					, col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
				);

			os << row;
		}	

		auto footage = 
			fmt::format( 
				"|={0:=<{1}}=|={0:=<{2}}=|={0:=<{3}}=|={0:=<{4}}=|={0:=<{5}}=|={0:=<{6}}=|\n" 

				, "", col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
			);

		os << footage;
	}

	// ### gnuplot 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot( GeometricPath const& gp )
	{
		std::string gnuplot_script_file_name( "dtg.gp" ); 
		// ##### open a file named {}
		std::ofstream gnuplot_script_file( gnuplot_script_file_name );

		// ##### write all lines of the triangulation in a datablock "line" and indexed by i
		std::for_each( edges(*this).first, edges(*this).second, 
			[&, i = 0, this] ( auto e ) mutable
			{
				auto const	s = source( e, *this);
				auto const	t = target( e, *this);
				auto const& ps = this->_points[ get( vertex_cell, *this, s )->source_index() ]; 
				auto const& pt = this->_points[ get( vertex_cell, *this, t )->source_index() ]; 

				gnuplot_script_file << 
					fmt::format( 
						"$line{0} << EOD \n"
						"{1} {2}\n" 
						"{3} {4}\n" 
						"EOD \n"
						, i++, ps.x(), ps.y(), pt.x(), pt.y()   
					);
			}	
		);

		// ##### write all points in a datablock named "points"
		gnuplot_script_file << 
			"$points << EOD \n";

		for (auto [it, it_end] = vertices(*this); it != it_end; ++it) {

			auto const& p =this->_points[ get( vertex_cell,*this, *it )->source_index() ]; 

			gnuplot_script_file << 
				fmt::format( 
					"{point_x} {point_y} {point_id} {pt_var} {ps_var} {lc_rgb_var} \n" 

					, "point_x"_a = p.x(), "point_y"_a = p.y(), "point_id"_a = *it 
					, "pt_var"_a = 7, "ps_var"_a = 1, "lc_rgb_var"_a = "0x66B2ff" 
				);
		}

		gnuplot_script_file << 
			"EOD \n";

		// ##### write all geometrical paths in a datablock named "path" and indexed by j
		std::for_each( std::begin( gp ), std::end( gp ), 
			[ &, j = 0, this ] ( auto const& path ) mutable
			{
				gnuplot_script_file 
					<< fmt::format( "$path{} << EOD \n", j++ );
				for ( auto const& point : path ) 
				{
					gnuplot_script_file << 
						fmt::format( 
							"{} {}\n" 
							, boost::geometry::get<0>(point), boost::geometry::get<1>(point)
						);
				}
				gnuplot_script_file 
					<< "EOD \n";
			}
		);
		// ##### write gnuplot commands and close the script file
		gnuplot_script_file << 
			"set title 'results' \n"
			"set style line 1 lt rgb 0x808080 \n"
			"set style line 2 lt rgb 0xB22222 lw 2 \n"
			"set style line 3 lt rgb 0x03C04A lw 1 \n"
			"plot 1/0 notitle \\\n";

		std::for_each( edges(*this).first, edges(*this).second, 
			[ &, i = 0, this ] ( auto e ) mutable
			{
				if ( get( boost::edge_index, *this, e ) == 3 ) { return; }	// do not plot all the nets
				gnuplot_script_file << 
					fmt::format( 
						", '$line{}' using 1:2 with lines ls {} title '' \\\n"
						, i++, ( get( boost::edge_index, *this, e ) == 1 ) ? 2 : 1 
					);
			}	
		);

		gnuplot_script_file <<
			", '$points' using 1:2:4:5:6 with points pt variable ps variable lc rgb variable title '' \\\n"
			", '$points' using 1:2:3 with labels offset char 0.5,0.5 title '' \\\n";

		for ( unsigned j = 0, j_end = gp.size(); j != j_end; ++j )
//		for ( unsigned j = 0, j_end = 9; j != j_end; ++j )
		{
			gnuplot_script_file <<
				fmt::format( 
					", '$path{}' using 1:2 with lines ls {} title '' \\\n"
					, j, 3 
				);
		}
		gnuplot_script_file.close();

		// ##### system call to gnuplot
		auto gnuplot_command = 
			fmt::format( "gnuplot -p {}", gnuplot_script_file_name );

		system( gnuplot_command.c_str() );

		// ##### system call to rm the gnuplot script 
		auto clean_commands =
			fmt::format( "rm  {}", gnuplot_script_file_name );

		system( clean_commands.c_str() );

	}

	// ### perturb 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	perturb( void )
	{
		std::default_random_engine 
			engine( std::chrono::system_clock::now().time_since_epoch().count() );
		std::uniform_int_distribution<std::size_t> 
			d( num_roots, num_leafs + num_roots - 1 );

        while (true)  
        {
            auto n = d( engine );
            _perturbs.push_back( n );
            if ( switch_mother( n ) ) { break; }
        }  
		//while ( !switch_mother( d( engine ) ) );
	}

	// ### switch_mother

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	bool planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	switch_mother( std::size_t u )
	{
		boost::filtered_graph											
			tri_edges( *this, detail::is_triangulation_edge( get( boost::edge_index, *this ) ) );

		auto M = get( vertex_mother, *this );
	
		auto const mu = get( vertex_mother, *this, u );
		auto edge_it_end = get( vertex_cell, *this, u )->incident_edge();
		while ( edge_it_end->twin()->cell()->source_index() != mu ) { edge_it_end = edge_it_end->prev(); }
		
		for ( auto edge_it = edge_it_end->prev(); edge_it != edge_it_end; edge_it = edge_it->prev() ) 
		{
			auto const x = edge_it->twin()->cell()->source_index();
			auto const [e1, b1] = edge( u, x, tri_edges );
			if ( !b1 ) continue;
			
			for ( auto v = x, w = get( vertex_mother, *this, v ); w != u; v = w, w = get( vertex_mother, *this, v ) )
			{
				if ( w == v )	// there is no loop
				{
					auto [e2, b2] = edge( u, mu, *this );
					put( boost::edge_index, *this, e2, 0 );
					put( boost::edge_index, *this, e1, 1 );
					put( vertex_mother, *this, u, x );
					return true;	
				}
			}
		} 
		return false;
	}

    // ### cost 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	double planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	cost( void )
	{
        unsigned cross_number = 0;
		boost::filtered_graph F( *this, detail::is_MST_edge( get( boost::edge_index, *this ) ) );
		boost::filtered_graph nets( *this, detail::is_net_edge( get( boost::edge_index, *this ) ) );
        
        /* version 1

		boost::depth_first_search( F, visitor( detail::accumulate_cross_number( num_roots, num_vertices( F ), nets, cross_number ) ) );

        */
        /* version 2

        // BFS traversal to obtain ordered foreset
        boost::adjacency_list<> 
            ordered_forest;
		for ( std::size_t i = 0; i != num_roots; ++i ) {
			boost::breadth_first_search( 
				F, i , 
	            visitor( detail::build_sketchable_forest( ordered_forest ) )
			);
		}
        // DFS traversal to obtain dest_order and compute cross_number
		boost::depth_first_search( 
			ordered_forest, 
			visitor( detail::compute_dest_ordering_and_cross_number( nets, num_vertices( ordered_forest ), num_roots, cross_number ) )
		);
        */

        // version 3

        // BFS traversal to obtain ordered foreset
        boost::adjacency_list<> 
            ordered_forest;
		for ( std::size_t i = 0; i != num_roots; ++i ) {
			boost::breadth_first_search( 
				F, i , 
	            visitor( detail::build_sketchable_forest( ordered_forest ) )
			);
		}
        // DFS traversal to obtain dest_order and compute cross_number
		boost::depth_first_search( 
			ordered_forest, 
			visitor( detail::compute_dest_ordering_and_cross_number2( nets, num_vertices( ordered_forest ), num_roots, cross_number ) )
		);
        //

		return cross_number;
	}


	// ## Private Methods of planar_orderd_forest_solution
	
	// ### construct_delaunay_trianguation
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	construct_delaunay_trianguation( void )	
	{
		auto v_cll_mp = boost::get( vertex_cell, *this );

		for (auto cell_it = _vd->cells().begin(); cell_it != _vd->cells().end(); ++cell_it) {

			cell_it->color( 1 );

			auto	u		= cell_it->source_index();
			auto	edge_it	= cell_it->incident_edge();

			do {

				if (auto adj_cell_it = edge_it->twin()->cell(); adj_cell_it->color() == 0) {

					auto v = adj_cell_it->source_index();

					add_edge( u, v, *this );
				}

				edge_it = edge_it->prev();

			} while (edge_it != cell_it->incident_edge());
		}
		for (auto cell_it = _vd->cells().begin(); cell_it != _vd->cells().end(); ++cell_it) {
			v_cll_mp[ cell_it->source_index() ] = cell_it;
		}
	}

	// ### mark_Planar_Ordered_Forest_edges
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	mark_Planar_Ordered_Forest_edges( void ) 
	{
		auto e_wgt_mp = boost::get( boost::edge_weight, *this ); 
		auto e_idx_mp = boost::get( boost::edge_index, *this ); 
		auto v_mom_mp = boost::get( vertex_mother, *this );
		auto v_deg_mp = boost::get( vertex_deg, *this );
		// #### set edge weight and find edges of MST by kruskal algorithm and mark them 1
		for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) {

			auto e = *it;
			auto u = source( e, *this );
			auto v = target( e, *this );
			if ( u > v ) std::swap( u, v );
			if ( v < num_roots && v - u == 1 ) 
			{ 
				e_wgt_mp[ e ] = 0; 
			}
			else 
			{
				e_wgt_mp[ e ] = boost::geometry::comparable_distance( _points[ u ], _points[ v ] );
			}
		}

		std::vector<
			typename boost::graph_traits<
				std::remove_reference_t<decltype( *this )>
			>::edge_descriptor 
		> 
			spanning_tree_edges;

		boost::kruskal_minimum_spanning_tree( *this, std::back_inserter( spanning_tree_edges ) );

		for (auto const& e: spanning_tree_edges) {
			e_idx_mp[ e ] = 1;								
		}

		// #### find boundary edges and mark them 2 
		for (std::size_t i = 0; i != num_roots; ++i) {
			auto [e, b] = edge( i, (i + 1) % num_roots, *this );
			boost::put( boost::edge_index, *this, e, 2 );
		}
		// #### set the vertex_mother_t property by BFS a filter graph of the MST
		boost::filtered_graph											// using std=c++17 feature. Cool!!!
			spannig_tree( *this, detail::is_MST_edge( e_idx_mp ) );

		for ( std::size_t i = 0; i != num_roots; ++i ) {
			v_mom_mp[ i ] = i;
			boost::breadth_first_search( 
				spannig_tree, i, 
				visitor( detail::set_mother_vertex_property( v_mom_mp ) )
			);
		}
		// #### find # of edges marked by 0 and set in the vertex_deg_t property in its adjacent edges
		for (auto [e_it, e_it_end] = edges( *this ); e_it != e_it_end; ++e_it) {

			if (e_idx_mp[ *e_it ] == 0) {

				++v_deg_mp[ source( *e_it, *this ) ];	
				++v_deg_mp[ target( *e_it, *this ) ];	

			}
		}
	}

	// ### mark_net_edges
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename NetContainer >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	mark_net_edges( NetContainer const& nets )
	{
		for ( auto const& n: nets ) {
			auto [e, b] = add_edge( std::get<0>( n ), std::get<1>( n ), *this );
			put( boost::edge_index, *this, e, 3 );
			put( boost::edge_weight, *this, e, std::numeric_limits<CoordinateType>::max() );
		}
	}
	

	// # circular_frame
	
	template < typename CircularFrameSolution >
	class circular_frame :
		public syc::topology::model::v_01::sketchable_forest< 
			typename CircularFrameSolution::coordinate_type, 
			typename CircularFrameSolution::point_type 
		>
	{
		// ## declaration

		private:
			using coordinate_type	= typename CircularFrameSolution::coordinate_type;
			using point_type		= typename CircularFrameSolution::point_type;
			using linestring_type	= typename CircularFrameSolution::linestring_type;
			using base				= syc::topology::model::v_01::sketchable_forest< point_type, coordinate_type >;

		// ## member data
		
		private:
			CircularFrameSolution&		_sol;

		// ## member functions 

		public:
			circular_frame( CircularFrameSolution& sol );
			auto the_basic_routing_algorithm( void );
			auto t_escape( void );
			auto t_escape2( void );
			auto realization( void );
            auto get_nets( void );

		private:
			bool clockwise_depth_first_route( typename base::topo_descriptor& h, std::size_t const u );
            
	};

	// ## Public Methods of circular_frame
	
	// ### constructor
	
	template < typename CircularFrameSolution >
	circular_frame< CircularFrameSolution >::
	circular_frame( CircularFrameSolution& sol ): _sol( sol )
	{
		boost::filtered_graph											
			spannig_tree( _sol, detail::is_MST_edge( get( boost::edge_index, _sol ) ) );
		for ( std::size_t i = 0; i != _sol.num_roots; ++i ) {
			boost::breadth_first_search( 
				spannig_tree, i , 
				visitor( detail::build_sketchable_forest( *this ) ) 
			);
		}
		this->make_sketch();
	}

	// ### the_basic_routing_algorithm
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	the_basic_routing_algorithm( void )
	{
		unsigned failure_cnt = 0; 

		boost::filtered_graph											
			nets( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );

		for ( auto [ ei, ei_end ] = edges( nets ); ei != ei_end; ++ei ) 
		{
			auto source_vertex = source( *ei, nets );					
			auto target_vertex = target( *ei, nets );					
			add_edge( source_vertex, target_vertex, this->nets );		// very ugly

			auto head = this->topo_instances( source_vertex )[ 0 ];	
			this->reset_slice_indexes();		// All the slices are not yet visited.
			if ( !clockwise_depth_first_route( head, target_vertex ) ) { ++failure_cnt; } 
		}	

		return failure_cnt;
	}

	// ### t_escape

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape( void )
	{
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::determine_net_priority_and_copmute_dest_ordering( *this, net, _sol.num_roots, dest_order ) ) 
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto T    = this->topo_instances( net_target );	
            auto input_it = std::find_if( std::begin( T ), std::end( T ), 
                [&]( auto t )
                {
                   return this->slice( h ) != this->slice( t ); 
                } 
            );

            this->reset_slice_indexes();            // All the slices are not yet visited.

            if ( input_it != std::end( T ) )
            {
                if ( !clockwise_depth_first_route( h, net_target ) ) { ++num_failed_nets; } 
            }
            else if ( dest_order[ net_source ][ 0 ] < 0 )       // should escape to the left of the current roots
            {
                auto anticlockwise_route_through = [&] ( auto& h, auto const& T )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; ++p != ori_mark; ) {
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                    for ( auto p = ori_mark; ++p != h; ) { 
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                };
                while ( !anticlockwise_route_through( h, T ) );
            }
            else
            {
                auto clockwise_route_through = [&] ( auto& h, auto const& T )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; --p != ori_mark; ) {
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                    for ( auto p = ori_mark; --p != h; ) { 
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                };
                while ( !clockwise_route_through( h, T ) );
            }
		}

		return num_failed_nets;
	}

	// ### t_escape2

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape2( void )
	{
        // version 1: shows a bug result from the make_slice operation positve edge CW route to next topoinstance
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::determine_net_priority_and_copmute_dest_ordering( *this, net, _sol.num_roots, dest_order ) )   
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto t    = this->topo_instances( net_target );	
            auto input_it = std::find_if( std::begin( t ), std::end( t ), 
                [&]( auto t )
                {
                   return this->slice( h ) != this->slice( t ); 
                } 
            );

            this->reset_slice_indexes();            // all the slices are not yet visited.

            if ( input_it != std::end( t ) )
            {
                if ( !clockwise_depth_first_route( h, net_target ) ) { ++num_failed_nets; } 
            }
            else if ( dest_order[ net_source ][ 0 ] < 0 )       // should escape to the left of the current roots
            {
                auto anticlockwise_route_through = [&] ( auto h, auto t )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; ++p != ori_mark; ) {
                        if ( p == t )
                        {
                            return this->make_slice( h, p );
                        } 
                        else if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) 
                            { 													
                                return this->make_slice( h, p );																	
                            }
                        }
                    }
                    for ( auto p = ori_mark; ++p != h; ) { 
                        if ( p == t )
                        {
                            return this->make_slice( h, p );
                        } 
                        else if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) 
                            { 													
                                return this->make_slice( h, p );																	
                            }
                        }
                    }
                };
                while ( h != t[0] )
                {
                    h = anticlockwise_route_through( h, t[ 0 ] );
                }
            }
            else
            {
                auto clockwise_route_through = [&] ( auto h, auto t )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; --p != ori_mark; ) {
                        if ( p == t )
                        {
                            return this->make_slice( h, p );
                        } 
                        else if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) { 													
                                e = this->make_slice( h, p );																	
                                return e;
                            }
                        }
                    }
                    for ( auto p = ori_mark; --p != h; ) { 
                        if ( p == t )
                        {
                            return this->make_slice( h, p );
                        } 
                        else if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) { 													
                                e = this->make_slice( h, p );																	
                                return e;
                            }
                        }
                    }
                };
                while ( h != t[0] )
                {
                    h = clockwise_route_through( h, t[0] );
                }
            }
		}

		return num_failed_nets;
    }

	// ### realization
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	realization( void ) 
	{
		auto R      = this->results();															    // every net should have a head before calling this method
		auto TCPs   = std::vector<std::vector<std::tuple<unsigned, unsigned, int, int>>>();	        // triangulation crossing paths
		for ( auto& path : R )  
		{						
			// step1: relax each of the paths
            auto dirty_patch = path.back();
			for ( bool relax = true; relax == true; ) {
				relax = false;
				auto it1 = path.begin();
				while ( true ) {
					auto it2 = std::next( it1 );
					if ( it2 == path.end() )	break;
					if ( it2->second == 0 )		break;
					if ( it2->second > 0 ) {
						auto cell_it = boost::get( vertex_cell, _sol, it1->first );
						auto vor_edge_it = cell_it->incident_edge();
						while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->next();	
						for ( ;    std::next( it2 )->first == vor_edge_it->next()->twin()->cell()->source_index()
								&& std::next( it2 )->second >= 0 && it1->second >=0 ;
							  vor_edge_it = vor_edge_it->next() ) {
							it2 = path.erase( it2 );
							relax = true;
						}
					}
					else if ( it2->second < 0 ) {
						auto cell_it = boost::get( vertex_cell, _sol, it1->first );
						auto vor_edge_it = cell_it->incident_edge();
						while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->prev();
						for ( ;    std::next( it2 )->first == vor_edge_it->prev()->twin()->cell()->source_index() 
								&& std::next( it2 )->second <= 0 && it1->second <=0 ;
							  vor_edge_it = vor_edge_it->prev() ) {
							it2 = path.erase( it2 );
							relax = true;
						}
					}
					it1 = it2;
				}
			}
            if ( path.back() != dirty_patch ) path.push_back( dirty_patch );
			// step2: realize each of the paths to a triangulation crossing paths
			std::vector<std::tuple<unsigned, unsigned, int, int>> X;
			if ( path.empty() )
			{
				TCPs.push_back( std::move( X ) );
			} 
			else
			{
				X.emplace_back( path.front().first, path.front().first, 0, 0 );
				auto	it1			= path.begin(), 
						it2			= std::next( it1 ),
						it3			= std::next( it2 );
				auto	vor_edge_it	= boost::get( vertex_cell, _sol, it1->first )->incident_edge();
				while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) 
				{
					vor_edge_it = vor_edge_it->next();
				}
				for	( ;it2->second != 0; ++it1, ++it2, ++it3 )
				{
					vor_edge_it = vor_edge_it->twin();
					if ( it2->second > 0 ) 
					{ 
						if ( it1->second < 0 )
						{
							X.emplace_back( it2->first, it1->first, it2->second, it1->second );
						}
						vor_edge_it = vor_edge_it->prev();
						for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
							  v != it3->first; 
							  vor_edge_it = vor_edge_it->prev(), v = vor_edge_it->twin()->cell()->source_index() )  
						{
							X.emplace_back( it2->first, v, it2->second, 0 );
						}
					} 
					else if ( it2->second < 0 ) 
					{
						if ( it1->second > 0 )
						{
							X.emplace_back( it2->first, it1->first, it2->second, it1->second );
						}
						vor_edge_it = vor_edge_it->next();
						for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
							  v != it3->first; 
							  vor_edge_it = vor_edge_it->next(), v = vor_edge_it->twin()->cell()->source_index() )  
						{
							X.emplace_back( it2->first, v, it2->second, 0 );
						}
					}
				}
				X.emplace_back( path.back().first, path.back().first, 0, 0 );
				TCPs.push_back( std::move( X ) );
			}
		}
        // step3: sort on edge
        auto GPds   = std::vector<std::vector<std::pair<typename std::list<int>::iterator,int>>>();
		for ( auto const& X : TCPs ) 
		{		
            auto gpd = std::vector<std::pair<typename std::list<int>::iterator,int>>();
			for ( auto const& t : X ) 
            { 
                auto [u, v, w, w2]  = t;
                auto n              = std::abs( w  );
                auto n2             = std::abs( w2 );
                if ( u == v ) { ; }
                else if ( u < v )    
                {
                    auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                    auto it = E[0].end(), it_end = it;
                    while ( ++it != it_end && *it < n ) { ; }
                    auto pos = E[0].insert( it, n );
                    gpd.emplace_back( pos, 0 );
                }
                else            
                {
                    if ( w2 == 0 )
                    {
                        auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                        auto it = E[1].end(), it_end = it;
                        while ( ++it != it_end && *it < n ) { ; }
                        auto pos = E[1].insert( it, n );
                        gpd.emplace_back( pos, 1 );
                    }
                    else
                    {
                        auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                        auto it = E[0].end(), it_end = it;
                        while ( ++it != it_end && *it < n2 ) { ; }
                        auto pos = E[0].insert( it, n2 );
                        gpd.emplace_back( pos, 0 );
                    }
                }
            }
            GPds.push_back( std::move( gpd ) );
        }
        // steps4: transform to geometrical paths
		auto GPs    = std::vector<linestring_type>();			            // geometrical paths
        for ( unsigned i = 0; i != TCPs.size(); ++i )
        {
			linestring_type gpath;
            for ( unsigned j = 0; j != TCPs[i].size(); ++j )
            {
                auto const& [u, v, w, w2] = TCPs[i][j];
				auto const& p0 = ( w2 != 0 && u > v )? _sol.point( v ): _sol.point( u );
				auto const& p1 = ( w2 != 0 && u > v )? _sol.point( u ): _sol.point( v );
				auto&       E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
				if ( u == v )
				{
				    auto const& p0 = _sol.point( u );
					auto const& x = boost::geometry::get<0>(p0);
					auto const& y = boost::geometry::get<1>(p0);
					boost::geometry::append( gpath, point_type( x, y ) );
				}
				else 
				{
                    double d = std::distance( E[GPds[i][j-1].second].begin(), GPds[i][j-1].first );
                    auto frac = (d + 1) / (E[0].size() + E[1].size() + 1);
					auto const& x = boost::geometry::get<0>(p0) - ( boost::geometry::get<0>(p0) - boost::geometry::get<0>(p1) ) * frac;
					auto const& y = boost::geometry::get<1>(p0) - ( boost::geometry::get<1>(p0) - boost::geometry::get<1>(p1) ) * frac;
					boost::geometry::append( gpath, point_type( x, y ) );
				}
			}
			GPs.push_back( std::move( gpath ) );
        }

		return GPs;
	}

	// ### get_nets
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	get_nets( void ) 
	{
        auto r = std::vector< std::array< std::size_t, 2 > >();
        for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei )
        {
            auto u = source( *ei, this->nets );
            auto v = target( *ei, this->nets );
            r.emplace_back( std::array< std::size_t, 2 >{ std::min( u , v ), std::max( u, v ) } );
        } 
        return r;
    }


	// ## Private Methods of circular_frame
	
	// ### clockwise_depth_first_route			
	
	template < typename CircularFrameSolution >
	bool circular_frame< CircularFrameSolution >::
	clockwise_depth_first_route( typename base::topo_descriptor& h, std::size_t const u ) 
	{
		// Clockwise seach for the first target vertex t of h that is within the same slice of h.
		// If sucessfully finds a t, route from h to t and move on to the next source vertex.
		auto m = this->slice( h ); 																		// Let m = the pointer to slice where h is in.
		m->index = 1; 																				// Set the color of slice m gray
		for ( auto const& t : this->topo_instances( u ) ) { 
			auto n = this->slice( t ); 																	// Let n = the pointer to slice where target pointed by t is in.
			if ( m == n ) {
				this->make_slice( h, t );	
				return true; 																		// Successfully find a topological routing sketch
			}
		}
		// EVENT POINT: cannot route from h to any of its targets in the current slice
		
		// Anticlockwise search for any candidate edges to go through
		std::stack<typename base::topo_descriptor> E;
		auto ori_mark = this->slice( h )->end();
		for ( auto p = h; ++p != ori_mark; ) {
			if ( this->attribution( p ) != 2 ) { 
				auto e = this->ordering_pair( p );															// Let e be the related edge of the edge pointed by p
				if ( 0 == this->slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
					E.push( p );																		// the edge pointed by pi is chosen as a candidate.
				}
			}
		}
		for ( auto p = ori_mark; ++p != h; ) { 
			if ( this->attribution( p ) != 2 ) { 
				auto e = this->ordering_pair( p );															// Let e be the related edge of the edge pointed by p
				if ( 0 == this->slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
					E.push( p );																		// the edge pointed by pi is chosen as a candidate.
				}
			}
		}
		while ( !E.empty() ) {
			auto inArrow = E.top();
			h = this->make_slice( h, inArrow );					
			if ( clockwise_depth_first_route( h, u ) ) { 
				return true; 
			}
			E.pop();
		}

		// EVENT POINT: In the current slice, there is no target to route to, and no edge go through. I.e Have to turn around, back to the previous slice.
		
		h = this->free_slice( h );

		return false;
	}


    // # package_data

    // ## detail of package_data

	namespace detail 
	{
        template < typename CoordinateType >
        struct io_pad
        {
            unsigned        id;
            CoordinateType  x;
            CoordinateType  y;
        };

        template < typename CoordinateType >
        struct net_type 
        {
            unsigned                                id;
            std::string                             name;
            std::vector< io_pad< CoordinateType > > pads;   // finger, bump, via
        };
    }

    // ## prototype of package_data

	template <
		typename CoordinateType = double, 
		typename PointType		= boost::polygon::point_data<CoordinateType>
	>
    class package_data  
    {
        // declartion

        using NetSeq = std::vector< detail::net_type< CoordinateType > >;
        
        // data members

        public:
            NetSeq              nets;
            CoordinateType      packagesizeX;
            CoordinateType      packagesizeY;
            CoordinateType      ballpitch; 
            CoordinateType      balldiameter;
            unsigned            balldimensionX;
            unsigned            balldimensionY; 
            CoordinateType      viadiamter; 
            CoordinateType      wirewidth; 
            CoordinateType      spacing;

        // constructor

        public:
            package_data( std::string const& ViaResults_FileName, std::string const& IoDrc_FileName );
        private:
            void read_via_results_file( std::string const& ViaResults_FileName );
            void read_io_drc_file( std::string const& IoDrc_FileName );
            void random_shuffle( detail::io_pad< CoordinateType >& p, std::uniform_real_distribution<CoordinateType>& u, std::default_random_engine& e );

        // gnuplot

        public:
            void gnuplot( std::string script_filename = "pckg.gp" );
        private:
            auto write_boudary_data( std::ofstream& script, char* color );
            auto write_bump_data( std::ofstream& script, char* color );
            auto write_finger_data( std::ofstream& script, char* color );
            auto write_via_data( std::ofstream& script, char* color );

        // generate_case

        public:
            auto generate_case( void );
            auto generate_fucked_up_case( void );
        private:

    };

    // ## constructor 

	template < typename CoordinateType, typename PointType >
    package_data< CoordinateType, PointType >:: 
    package_data( std::string const& ViaResults_FileName, std::string const& IoDrc_FileName )
    {
        read_via_results_file( ViaResults_FileName );
        read_io_drc_file( IoDrc_FileName );
    }

    // ### read_via_results_file

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    read_via_results_file( std::string const& ViaResults_FileName )
    {
        auto infile     = std::ifstream( ViaResults_FileName );
        auto token      = std::string();
        auto c          = char(); 
        auto u          = std::uniform_real_distribution<CoordinateType>( -10, 10 );
        auto e          = std::default_random_engine( 0 );

        while ( getline( infile, token ) )
        {
            auto reg    = std::regex( "([0-9]+)" );
            auto sm     = std::smatch();
            if ( !std::regex_match( token, sm, reg ) ) continue;

            auto i = std::stoul( sm[1] ); 

            auto t = detail::net_type< CoordinateType >(); 
            t.id = i;

            auto read = [&] ( auto keyword )
            {
                assert( infile >> token );

                if ( token == keyword && keyword == "Net" )
                {
                    assert( infile >> c >> t.name );
                }
                else if ( token == keyword && keyword == "Finger" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.id >> c >> p.x >> p.y >> c );
                    t.pads.push_back( std::move( p ) );
                }
                else if ( token == keyword && keyword == "Bump" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.id >> c >> p.x >> p.y >> c );
                    random_shuffle( p, u, e );
                    t.pads.push_back( std::move( p ) );
                }
                else if ( token == keyword && keyword == "Via" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.x >> p.y >> c );
                    p.id = i;
                    t.pads.push_back( std::move( p ) );
                }
                else
                {
                    assert( false );
                }                    
            };

            read( "Net"     );
            read( "Finger"  );
            read( "Bump"    );
            read( "Via"     );

            nets.push_back( std::move( t ) );
        }
    }

    // ### read_io_drc_file

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    read_io_drc_file( std::string const& IoDrc_FileName )
    {
        auto infile     = std::ifstream( IoDrc_FileName );
        auto token      = std::string();

        infile >> token >> packagesizeX;
        infile >> token >> packagesizeY;
        infile >> token >> ballpitch;
        infile >> token >> balldiameter;
        infile >> token >> balldimensionX;
        infile >> token >> balldimensionY;
        infile >> token >> viadiamter; 
        infile >> token >> wirewidth; 
        infile >> token >> spacing; 
    }

    // ### random_shuffle
	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    random_shuffle( detail::io_pad< CoordinateType >& p, std::uniform_real_distribution<CoordinateType>& u, std::default_random_engine& e )
    {
        p.x += u( e );
        p.y += u( e );
    }


    // ## gnuplot 

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    gnuplot( std::string script_filename )
    {
        auto script            = std::ofstream( script_filename );

        auto boundary_plot_cmd = write_boudary_data( script, "0x808080" );  // color palette1: https://www.pinterest.com/pin/25755029108705685/?d=t&mt=login
        auto bump_plot_cmd     = write_bump_data( script, "0xB7C5C8" );     // color palette2: https://www.pinterest.com/pin/222083825344676565/
        auto vias_plot_cmd     = write_via_data( script, "0xE49969" );
        auto fingers_plot_cmd  = write_finger_data( script, "0xC9C27F" );

        script  << fmt::format( "plot {}    \\\n", boundary_plot_cmd    )
                << fmt::format( ",    {}    \\\n", bump_plot_cmd        )
                << fmt::format( ",    {}    \\\n", vias_plot_cmd        )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[0]  )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[1]  )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[2]  );

        script.close();

        system( fmt::format( "gnuplot -p {}", script_filename ).c_str() ); 
        system( fmt::format( "rm {}"        , script_filename ).c_str() ); 
    }

    // ### write_boudary

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_boudary_data( std::ofstream& script, char* color )
    {
        auto x = packagesizeX / 2;
        auto y = packagesizeY / 2;
		script << "set style line 1 lt rgb " << color << "\n";
		script << "$boundary << EOD \n";
		script << fmt::format( 
            "{x1} {y1} \n" 
            "{x2} {y2} \n" 
            "{x3} {y3} \n" 
            "{x4} {y4} \n" 
            "{x1} {y1} \n" 
            , fmt::arg( "x1",  x ), fmt::arg( "y1",  y )
            , fmt::arg( "x2", -x ), fmt::arg( "y2",  y )
            , fmt::arg( "x3", -x ), fmt::arg( "y3", -y )
            , fmt::arg( "x4",  x ), fmt::arg( "y4", -y )
        ); 
		script << "EOD \n";
        return "'$boundary' using 1:2 with lines ls 1 notitle";        
    }

    // ### write_bump_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_bump_data( std::ofstream& script, char* color )
    {
        auto x = ballpitch * (balldimensionX - 1) / 2;
        auto y = ballpitch * (balldimensionY - 1) / 2;
        auto r = balldiameter / 2;
		script << "$bumps << EOD \n";
        for ( unsigned i = 0; i != balldimensionX; ++i )
        {
            for ( unsigned j = 0; j != balldimensionY; ++j )
            {
                script << fmt::format( 
                    "{x} {y} {radius} {color} \n" 
                    , fmt::arg( "x"             , x - i * ballpitch )
                    , fmt::arg( "y"             , y - j * ballpitch )
                    , fmt::arg( "radius"        , r                 )
                    , fmt::arg( "color"         , color             )
                );
            }
        }
		script << "EOD \n";
        return "'$bumps' using 1:2:3:4 with circles lc rgb var notitle";
    }

    // ### write_via_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_via_data( std::ofstream& script, char* color )
    {
		script << "$vias << EOD \n";
        std::for_each( nets.begin(), nets.end(),    
            [ & ] ( auto n ) mutable
            {
                script << fmt::format( 
                    "{x} {y} {pt_var} {ps_var} {lc_rgb_var} \n" 
                    , fmt::arg( "x"             , n.pads[2].x   )
                    , fmt::arg( "y"             , n.pads[2].y   )
                    , fmt::arg( "pt_var"        , 7             )
                    , fmt::arg( "ps_var"        , 1             )
                    , fmt::arg( "lc_rgb_var"    , color         )
                );
            }
        );
		script << "EOD \n";
		return "'$vias' using 1:2:3:4:5 with points pt variable ps variable lc rgb variable notitle";
    }

    // ### write_finger_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_finger_data( std::ofstream& script, char* color )
    {
        auto c1 = "0xB7C5C8";
		script << "set style line 2 lt rgb " << color << "\n";
		script << "$fingers << EOD \n";
        for ( auto const& n : nets )
        {
            script << fmt::format( 
                "{x} {y} {id} {pt_var} {ps_var} {lc_rgb_var} \n" 
                , fmt::arg( "x"             , n.pads[0].x   )
                , fmt::arg( "y"             , n.pads[0].y   )
                , fmt::arg( "id"            , n.pads[0].id  )
                , fmt::arg( "pt_var"        , 7             )
                , fmt::arg( "ps_var"        , 1             )
                , fmt::arg( "lc_rgb_var"    , color         )
            );
        }
		script << "EOD \n";
		return std::array< std::string, 3 > 
        {   "'$fingers' using 1:2:4:5:6 with points pt variable ps variable lc rgb variable notitle",	
			"'$fingers' using 1:2:3 with labels offset char 0.5,0.5 notitle",
            "'$fingers' using 1:2 with lines ls 2 notitle"
        };
    }


    // ## generate_case

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_case( void )
    {
        auto _roots         = std::vector< PointType >();
        auto _inner_bumps   = std::vector< PointType >();
        auto _outer_bumps   = std::vector< PointType >();
        auto _inner_nets    = std::vector< std::array< unsigned, 2 > >(); 
        auto _outer_nets    = std::vector< std::array< unsigned, 2 > >(); 

        for ( auto const& e : nets )
        {
            _roots.emplace_back( e.pads[0].x, e.pads[0].y );
        }
        _roots.push_back( _roots.front() );

        for ( unsigned i = 0; i != nets.size(); ++i )
        {
            auto p = PointType( nets[ i ].pads[ 1 ].x, nets[ i ].pads[ 1 ].y );
            if ( boost::geometry::within( p, _roots ) )
            {
                _inner_nets.emplace_back( 
                    std::array{ (unsigned)(nets.size() + _inner_bumps.size()), i } 
                );
                _inner_bumps.push_back( std::move( p ) );
            }
            else
            {
                _outer_nets.emplace_back( 
                    std::array{ (unsigned)(nets.size() + _outer_bumps.size()), i } 
                );
                _outer_bumps.push_back( std::move( p ) );
            }
        }

        _roots.pop_back();

        return std::tuple( _roots, _inner_bumps, _inner_nets, _outer_bumps, _outer_nets );
    }


}
namespace syc::topology::v_03 
{
	using namespace fmt::literals;

    // # package_data

    // ## detail of package_data

	namespace detail 
	{
        // === io_pad ===
        template < typename CoordinateType >
        struct io_pad
        {
            unsigned        id;
            CoordinateType  x;
            CoordinateType  y;
        };

        // === net_type ===
        template < typename CoordinateType >
        struct net_type 
        {
            unsigned                                id;
            std::string                             name;
            std::vector< io_pad< CoordinateType > > pads;   // finger, bump, via
        };
    }

    // ## prototype of package_data

	template <
		typename CoordinateType = double, 
		typename PointType		= boost::polygon::point_data<CoordinateType>
	>
    class package_data  
    {
        // declartion

        using coordinate_type   = CoordinateType;
        using point_type        = PointType;
        using NetSeq            = std::vector< detail::net_type< CoordinateType > >;
        
        // data members

        public:
            NetSeq              nets;
            CoordinateType      packagesizeX;
            CoordinateType      packagesizeY;
            CoordinateType      ballpitch; 
            CoordinateType      balldiameter;
            unsigned            balldimensionX;
            unsigned            balldimensionY; 
            CoordinateType      viadiamter; 
            CoordinateType      wirewidth; 
            CoordinateType      spacing;

        // constructor

        public:
            package_data( std::string const& ViaResults_FileName, std::string const& IoDrc_FileName );
        private:
            void read_via_results_file( std::string const& ViaResults_FileName );
            void read_io_drc_file( std::string const& IoDrc_FileName );
            void random_shuffle( detail::io_pad< CoordinateType >& p, std::uniform_real_distribution<CoordinateType>& u, std::default_random_engine& e );
            void reorder_nets_by_computing_cross_product_of_finger_position( void );

        // gnuplot

        public:
            void gnuplot( std::string script_filename = "pckg.gp" );
        private:
            auto write_boudary_data( std::ofstream& script, char* color );
            auto write_bump_data( std::ofstream& script, char* color );
            auto write_finger_data( std::ofstream& script, char* color );
            auto write_via_data( std::ofstream& script, char* color );

        // generate_case (using_bumps or using vias)

        public:
            auto generate_case_using_bumps( void );
            auto generate_case_using_vias( void ); 
        private:
            auto generate_case( int );

        // generate a fucked up case to graguate

        public:
            auto generate_fucked_up_case_using_bumps( void );
            auto generate_fucked_up_case_using_vias( void ); 
        private:
            auto generate_fucked_up_case( int BorV );
	        template < typename PairPPointNetindex >
            auto debugging_in_GIF( PairPPointNetindex const& Proj_P, std::string GIF_name );
            auto cal_c_P( void );
            template < typename MultiPoint >
            auto proj_to_bdry( PointType const& c, MultiPoint const& P );
    
    };

    // ## constructor 

	template < typename CoordinateType, typename PointType >
    package_data< CoordinateType, PointType >:: 
    package_data( std::string const& ViaResults_FileName, std::string const& IoDrc_FileName )
    {
        read_via_results_file( ViaResults_FileName );
        read_io_drc_file( IoDrc_FileName );
        //reorder_nets_by_computing_cross_product_of_finger_position();   // turn off when running 99io
    }

    // ### read_via_results_file

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    read_via_results_file( std::string const& ViaResults_FileName )
    {
        auto infile     = std::ifstream( ViaResults_FileName );
        auto token      = std::string();
        auto c          = char(); 
        auto u          = std::uniform_real_distribution<CoordinateType>( -10, 10 );
        auto e          = std::default_random_engine( 0 );

        while ( getline( infile, token ) )
        {
            auto reg    = std::regex( "([0-9]+)" );
            auto sm     = std::smatch();
            if ( !std::regex_match( token, sm, reg ) ) continue;

            auto i = std::stoul( sm[1] ); 

            auto t = detail::net_type< CoordinateType >(); 
            t.id = i;

            auto read = [&] ( auto keyword )
            {
                assert( infile >> token );

                if ( token == keyword && keyword == "Net" )
                {
                    assert( infile >> c >> t.name );
                }
                else if ( token == keyword && keyword == "Finger" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.id >> c >> p.x >> p.y >> c );
                    t.pads.push_back( std::move( p ) );
                }
                else if ( token == keyword && keyword == "Bump" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.id >> c >> p.x >> p.y >> c );
                    random_shuffle( p, u, e );
                    t.pads.push_back( std::move( p ) );
                }
                else if ( token == keyword && keyword == "Via" )
                {
                    auto p = detail::io_pad< CoordinateType >();
                    assert( infile >> c >> p.x >> p.y >> c );
                    p.id = i;
                    t.pads.push_back( std::move( p ) );
                }
                else
                {
                    assert( false );
                }                    
            };

            read( "Net"     );
            read( "Finger"  );
            read( "Bump"    );
            read( "Via"     );

            nets.push_back( std::move( t ) );
        }
    }

    // ### read_io_drc_file

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    read_io_drc_file( std::string const& IoDrc_FileName )
    {
        auto infile     = std::ifstream( IoDrc_FileName );
        auto token      = std::string();

        infile >> token >> packagesizeX;
        infile >> token >> packagesizeY;
        infile >> token >> ballpitch;
        infile >> token >> balldiameter;
        infile >> token >> balldimensionX;
        infile >> token >> balldimensionY;
        infile >> token >> viadiamter; 
        infile >> token >> wirewidth; 
        infile >> token >> spacing; 
    }

    // ### random_shuffle

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    random_shuffle( detail::io_pad< CoordinateType >& p, std::uniform_real_distribution<CoordinateType>& u, std::default_random_engine& e )
    {
        p.x += u( e );
        p.y += u( e );
    }

    // ### reorder_nets_by_computing_cross_product_of_finger_position

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    reorder_nets_by_computing_cross_product_of_finger_position( void )
    {
        std::sort( nets.begin(), nets.end(), 
            [&] ( auto const& n1, auto const& n2 )
            {
                auto p1 = point_type( n1.pads[0].x, n1.pads[0].y );
                auto p2 = point_type( n2.pads[0].x, n2.pads[0].y );
                auto p3 = boost::geometry::cross_product( p1, p2 );
                return boost::geometry::get<0>( p3 ) > 0;
            }
        );
    }


    // ## gnuplot 

	template < typename CoordinateType, typename PointType >
    void package_data< CoordinateType, PointType >:: 
    gnuplot( std::string script_filename )
    {
        auto script            = std::ofstream( script_filename );

        auto boundary_plot_cmd = write_boudary_data( script, "0x808080" );  // color palette1: https://www.pinterest.com/pin/25755029108705685/?d=t&mt=login
        auto bump_plot_cmd     = write_bump_data( script, "0xB7C5C8" );     // color palette2: https://www.pinterest.com/pin/222083825344676565/
        auto vias_plot_cmd     = write_via_data( script, "0xE49969" );
        auto fingers_plot_cmd  = write_finger_data( script, "0xC9C27F" );

        script  << fmt::format( "plot {}    \\\n", boundary_plot_cmd    )
                << fmt::format( ",    {}    \\\n", bump_plot_cmd        )
                << fmt::format( ",    {}    \\\n", vias_plot_cmd        )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[0]  )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[1]  )
                << fmt::format( ",    {}    \\\n", fingers_plot_cmd[2]  );

        script.close();

        system( fmt::format( "gnuplot -p {}", script_filename ).c_str() ); 
        system( fmt::format( "rm {}"        , script_filename ).c_str() ); 
    }

    // ### write_boudary

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_boudary_data( std::ofstream& script, char* color )
    {
        auto x = packagesizeX / 2;
        auto y = packagesizeY / 2;
		script << "set style line 1 lt rgb " << color << "\n";
		script << "$boundary << EOD \n";
		script << fmt::format( 
            "{x1} {y1} \n" 
            "{x2} {y2} \n" 
            "{x3} {y3} \n" 
            "{x4} {y4} \n" 
            "{x1} {y1} \n" 
            , fmt::arg( "x1",  x ), fmt::arg( "y1",  y )
            , fmt::arg( "x2", -x ), fmt::arg( "y2",  y )
            , fmt::arg( "x3", -x ), fmt::arg( "y3", -y )
            , fmt::arg( "x4",  x ), fmt::arg( "y4", -y )
        ); 
		script << "EOD \n";
        return "'$boundary' using 1:2 with lines ls 1 notitle";        
    }

    // ### write_bump_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_bump_data( std::ofstream& script, char* color )
    {
        auto x = ballpitch * (balldimensionX - 1) / 2;
        auto y = ballpitch * (balldimensionY - 1) / 2;
        auto r = balldiameter / 2;
		script << "$bumps << EOD \n";
        for ( unsigned i = 0; i != balldimensionX; ++i )
        {
            for ( unsigned j = 0; j != balldimensionY; ++j )
            {
                script << fmt::format( 
                    "{x} {y} {radius} {color} \n" 
                    , fmt::arg( "x"             , x - i * ballpitch )
                    , fmt::arg( "y"             , y - j * ballpitch )
                    , fmt::arg( "radius"        , r                 )
                    , fmt::arg( "color"         , color             )
                );
            }
        }
		script << "EOD \n";
        return "'$bumps' using 1:2:3:4 with circles lc rgb var notitle";
    }

    // ### write_via_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_via_data( std::ofstream& script, char* color )
    {
		script << "$vias << EOD \n";
        std::for_each( nets.begin(), nets.end(),    
            [ & ] ( auto n ) mutable
            {
                script << fmt::format( 
                    "{x} {y} {pt_var} {ps_var} {lc_rgb_var} \n" 
                    , fmt::arg( "x"             , n.pads[2].x   )
                    , fmt::arg( "y"             , n.pads[2].y   )
                    , fmt::arg( "pt_var"        , 7             )
                    , fmt::arg( "ps_var"        , 1             )
                    , fmt::arg( "lc_rgb_var"    , color         )
                );
            }
        );
		script << "EOD \n";
		return "'$vias' using 1:2:3:4:5 with points pt variable ps variable lc rgb variable notitle";
    }

    // ### write_finger_data

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    write_finger_data( std::ofstream& script, char* color )
    {
        auto c1 = "0xB7C5C8";
		script << "set style line 2 lt rgb " << color << "\n";
		script << "$fingers << EOD \n";
        for ( auto const& n : nets )
        {
            script << fmt::format( 
                "{x} {y} {id} {pt_var} {ps_var} {lc_rgb_var} \n" 
                , fmt::arg( "x"             , n.pads[0].x   )
                , fmt::arg( "y"             , n.pads[0].y   )
                , fmt::arg( "id"            , n.pads[0].id  )
                , fmt::arg( "pt_var"        , 7             )
                , fmt::arg( "ps_var"        , 1             )
                , fmt::arg( "lc_rgb_var"    , color         )
            );
        }
		script << "EOD \n";
		return std::array< std::string, 3 > 
        {   "'$fingers' using 1:2:4:5:6 with points pt variable ps variable lc rgb variable notitle",	
			"'$fingers' using 1:2:3 with labels offset char 0.5,0.5 notitle",
            "'$fingers' using 1:2 with lines ls 2 notitle"
        };
    }


    // ## generate_case (using_bumps or using vias)

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_case_using_bumps( void )
    {
        return generate_case( 1 );
    }
	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_case_using_vias( void )
    {
        return generate_case( 2 );
    }
	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_case( int BorV )
    {
        auto _roots         = std::vector< PointType >();
        auto _inner_points  = std::vector< PointType >();
        auto _outer_points  = std::vector< PointType >();
        auto _inner_nets    = std::vector< std::array< unsigned, 2 > >(); 
        auto _outer_nets    = std::vector< std::array< unsigned, 2 > >(); 

        for ( auto const& e : nets )
        {
            _roots.emplace_back( e.pads[0].x, e.pads[0].y );
        }

        _roots.push_back( _roots.front() );

        assert( boost::geometry::is_simple( _roots ) );

        for ( unsigned i = 0; i != nets.size(); ++i )
        {
            auto p = PointType( nets[ i ].pads[ BorV ].x, nets[ i ].pads[ BorV ].y );
            if ( boost::geometry::within( p, _roots ) )
           {
                _inner_nets.emplace_back( 
                    std::array{ (unsigned)(nets.size() + _inner_points.size()), i } 
                );
                _inner_points.push_back( std::move( p ) );
            }
            else
            {
                _outer_nets.emplace_back( 
                    std::array{ (unsigned)(nets.size() + _outer_points.size()), i } 
                );
                _outer_points.push_back( std::move( p ) );
            }
        }

        _roots.pop_back();

        return std::tuple( _roots, _inner_points, _inner_nets, _outer_points, _outer_nets );

    }

	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_fucked_up_case_using_bumps( void )
    {
        return generate_fucked_up_case( 1 );
    }
	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_fucked_up_case_using_vias( void )
    {
        return generate_fucked_up_case( 2 );
    }
	template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    generate_fucked_up_case( int BorV )
    {
        // return c: the centroid of the fingers
        //        P: mathematic vectors from the centroid to fingers indexed by net
        auto [c, P] = cal_c_P(); 

        // first:   projected finger position
        // second:  its net index
        // boundary points included.
        auto Proj_P = proj_to_bdry( c, P );

        // sorted Proj_P by theta
        std::sort( Proj_P.begin(), Proj_P.end(), 
            []( auto const p1, auto const p2 )
            {
                auto theta1 = std::atan2( boost::geometry::get<1>( p1.first ), boost::geometry::get<0>( p1.first ) ); 
                auto theta2 = std::atan2( boost::geometry::get<1>( p2.first ), boost::geometry::get<0>( p2.first ) ); 
                return  theta1 < theta2;     
            } 
        );

//      debugging_in_GIF( Proj_P, "sorted_finger.gif" );

        auto _roots         = std::vector< PointType >();
        auto _inner_points  = std::vector< PointType >();
        auto _inner_nets    = std::vector< std::array< unsigned, 2 > >(); 

        for ( unsigned i = 0; i != Proj_P.size(); ++i )
        {
            _roots.push_back( Proj_P[i].first );
            auto j = Proj_P[i].second;
            if ( j < P.size() ) 
            {
                _inner_nets.emplace_back( 
                    std::array{ (unsigned)(Proj_P.size() + _inner_points.size()), i } 
                );
                _inner_points.push_back( PointType( nets[ j ].pads[ BorV ].x, nets[ j ].pads[ BorV ].y ) );
            }
        }            

        return  std::tuple( _roots, _inner_points, _inner_nets );
    }
    template < typename CoordinateType, typename PointType >
    auto package_data< CoordinateType, PointType >:: 
    cal_c_P( void )
    {
        auto P = boost::geometry::model::multi_point<PointType>();
        for ( auto const& e : nets )
        {
            boost::geometry::append( P, PointType( e.pads[0].x, e.pads[0].y ) );
        }

        auto c = PointType();
        boost::geometry::centroid( P, c );

        for ( auto& p : P )
        {
            boost::geometry::subtract_point( p, c );    // p = p - c;
        }
        return std::pair( c, P );
    }
	template < typename CoordinateType, typename PointType >
	template < typename MultiPoint >
    auto package_data< CoordinateType, PointType >:: 
    proj_to_bdry( PointType const& c, MultiPoint const& P )
    {
        auto ret = std::vector< std::pair< PointType, unsigned > >();

        using SegmentType = boost::geometry::model::segment< PointType >;
        auto B = std::vector< SegmentType >();
        auto x0 = (balldimensionX) * ballpitch / 2;
        auto y0 = (balldimensionY) * ballpitch / 2;
        B.push_back( SegmentType{ {  x0,  y0 }, { -x0,  y0 } } );
        B.push_back( SegmentType{ { -x0,  y0 }, { -x0, -y0 } } );
        B.push_back( SegmentType{ { -x0, -y0 }, {  x0, -y0 } } );
        B.push_back( SegmentType{ {  x0, -y0 }, {  x0,  y0 } } );

        
        for ( unsigned i = 0, i_end = P.size(); i != i_end; ++i )
        {
            auto output = std::vector< PointType >();
            auto k = c;
            boost::geometry::add_point( k, P[i] );
            for ( bool flag = true; flag; )
            {
                boost::geometry::add_point( k, P[i] );
                auto CK = boost::geometry::model::segment< PointType >{ c, k };

                for ( unsigned j = 0; j != 4; ++j )
                {
                    boost::geometry::intersection( CK, B[j], output );  
                    if ( !output.empty() )
                    {
                        ret.emplace_back( output.front(), i );
                        flag = false;
                        break;
                    }
                }
            }
        }

        ret.push_back( std::pair< PointType, unsigned >( PointType( x0,  y0), P.size() + 0 ) );
        ret.push_back( std::pair< PointType, unsigned >( PointType(-x0,  y0), P.size() + 1 ) );
        ret.push_back( std::pair< PointType, unsigned >( PointType(-x0, -y0), P.size() + 2 ) );
        ret.push_back( std::pair< PointType, unsigned >( PointType( x0, -y0), P.size() + 3 ) );

        return ret;
    }
	
	template < typename CoordinateType, typename PointType >
	template < typename PairPPointNetindex >
    auto package_data< CoordinateType, PointType >:: 
    debugging_in_GIF( PairPPointNetindex const& Proj_P, std::string GIF_name )
    {
        auto data_file_name   = "sorted_finger.txt";
        auto plot_script_name = "sorted_finger.gp";

        auto* outfile = fopen( data_file_name, "w" );
        {
            for ( unsigned i = 0; i != Proj_P.size(); ++i )
            {
                fprintf( outfile, "%lf %lf \n", 0, 0 );
                fprintf( outfile, "%lf %lf \n\n", boost::geometry::get<0>( Proj_P[i].first ), boost::geometry::get<1>( Proj_P[i].first ) );
            }
        }
        fclose( outfile );
        auto* gnuplot = fopen( plot_script_name, "w" );
        {
            auto const x = (balldimensionX) * ballpitch / 2;
            auto const y = (balldimensionY) * ballpitch / 2;
            fmt::print( gnuplot, "set term gif animate\n"                                              );
            fmt::print( gnuplot, "set output '{}'\n", GIF_name                                         );
            fmt::print( gnuplot, "set xrange [-{}:{}]\n", packagesizeX/2, packagesizeY/2               );
            fmt::print( gnuplot, "set yrange [-{}:{}]\n", packagesizeX/2, packagesizeY/2               );
            fmt::print( gnuplot, "$inner_boundary << EOD\n"                                            );
            fmt::print( gnuplot, "{} {}\n",  x ,  y                                                    );  
            fmt::print( gnuplot, "{} {}\n", -x ,  y                                                    );  
            fmt::print( gnuplot, "{} {}\n", -x , -y                                                    );  
            fmt::print( gnuplot, "{} {}\n",  x , -y                                                    );  
            fmt::print( gnuplot, "EOD\n"                                                               );
            fmt::print( gnuplot, "do for [i=0:175]{{\n"                                                );
            fmt::print( gnuplot, "    plot '{}' every :::0::i with lines notitle \\\n", data_file_name );
            fmt::print( gnuplot, "    , $inner_boundary with points notitle        \n"                 );
            fmt::print( gnuplot, "}}\n"                                                                );
            fmt::print( gnuplot, "set output\n"                                                        );
        }
        fclose( gnuplot );

        system( fmt::format( "gnuplot '{}'" , plot_script_name       ).c_str() );
        system( fmt::format( "rm {}"        , data_file_name         ).c_str() );
        system( fmt::format( "rm {}"        , plot_script_name       ).c_str() );
        system( fmt::format( "mv {} {}"     , GIF_name, "~/Desktop/" ).c_str() );
    }


	// # planar_orderd_forest_solution

    // ## detail        of planar_orderd_forest_solution
    
    namespace detail 
    { 
		// ### predicate: is_triangulation_edge
		
		template <typename EdgeIndexMap>
		struct is_triangulation_edge 
		{
			public:
				is_triangulation_edge() = default;
				is_triangulation_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 0;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ### predicate: is_MST_edge
		
		template <typename EdgeIndexMap>
		struct is_MST_edge 
		{
			public:
				is_MST_edge() = default;
				is_MST_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 1;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ### predicate: is_boundary_edge
		
		template <typename EdgeIndexMap>
		struct is_boundary_edge 
		{
			public:
				is_boundary_edge() = default;
				is_boundary_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 2;
				}
			private:
				EdgeIndexMap I;
		};
		
		// ### predicate: is_net_edge
		
		template <typename EdgeIndexMap>
		struct is_net_edge 
		{
			public:
				is_net_edge() = default;
				is_net_edge( EdgeIndexMap m ): I( m ) {}
				template <typename Edge>
				bool operator()( const Edge& e ) const 
				{
					return I[ e ] == 3;
				}
			private:
				EdgeIndexMap I;
		};
		
        // ### compute_dest_ordering_and_cross_number2

		template < typename Net >
        struct compute_dest_ordering_and_cross_number2 : public boost::default_dfs_visitor 
        {
            // declarition
            using DestOrder = std::vector< std::pair< int, std::size_t > >;
            using SortingStack = std::vector< std::list< std::size_t > >;
            // data member
            Net const&              N;          // IN:      nets
            int                     B;          // IN:      number of roots on the boudary
            int                     R;          // UTIL:    current root
            DestOrder               D;          // UTIL:    dest_order
            SortingStack            S;          // UTIL:    sort to compute cross number
            unsigned&               cross_n;    // OUT:     cross_number
            // constructor
            compute_dest_ordering_and_cross_number2( Net const& n, int num_vertices, int num_roots, unsigned& cn ):
                N( n ), B( num_roots ), D( num_vertices ), cross_n( cn ) {}
            // start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
            // finish_vertex
			template < typename Vertex, typename Graph >
			void finish_vertex( Vertex u, Graph const& g ) 
			{
                compute_D( u );

                if ( out_degree( u, g ) == 0  ) {
                    if ( u >= B ) { S.emplace_back( 1, u ); }   // it is a leaf
                }
                else
                {
                    auto comp = [&]( auto f2, auto f1 ) 
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        return d2 < d1;
                    };
                    unsigned f1_left;
                    auto comp2 = [&]( auto f2, auto f1 )
                    {
                        auto d2 = std::get<0>( D[f2] ) - long( std::get<1>( D[f2] ) );
                        auto d1 = std::get<0>( D[f1] ) - long( std::get<1>( D[f1] ) );
                        if ( d2 < d1 )
                        {
                            cross_n += f1_left;
                            return true;
                        }
                        else
                        {
                            --f1_left;
                            return false;
                        }
                    };
                    std::list<std::size_t> L( std::move( S.back() ) );
                    S.pop_back(); 
                    auto merge_and_compute_cross_n = [&] ( void )
                    {
                        f1_left = L.size();
                        L.merge( std::move( S.back() ), comp2 );
                        S.pop_back(); 
                    };
                    for ( unsigned i = 0, i_end = out_degree( u, g ) - 1; i != i_end; ++i ) 
                    {
                        merge_and_compute_cross_n();
                    }
                    if ( u < B ) 
                    {
                        auto l = std::get<0>( D[ L.front() ] ) - long( std::get<1>( D[ L.front() ] ) );
                        auto r = std::get<0>( D[ L.back() ] ) - long( std::get<1>( D[ L.back() ] ) );
                        auto sz = L.size();
                        auto m = r - l + 1 - sz;
                        cross_n += m * (sz + 1) / 2;
                    }
                    else 
                    {
                        L.merge( std::list<std::size_t>( 1, u ), comp );
                    }
                    S.push_back( std::move( L ) );
                }
            }
            // compute_D
			template < typename Vertex >
            void compute_D( Vertex u )
            {
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    std::get<0>( D[ u ] ) = ( Dcw < Dccw )? Dcw : -Dccw;
                    std::get<1>( D[ u ] ) = R;
                }
                else
                {
                    std::get<0>( D[ u ] ) = 0;
                    std::get<1>( D[ u ] ) = R;
                } 
            }
        };

		// ### build_sketchable_forest
		
		template < typename sketable_forest >
		class build_sketchable_forest: public boost::default_bfs_visitor 
		{
			// constructor 
			public:
				build_sketchable_forest( sketable_forest& sf ): _SF( sf ) {}
			// data member
			private:
				sketable_forest& _SF;
			// events
			public:	
				// examine_vertex
				template <typename Vertex, typename Graph>	
					void examine_vertex( Vertex u, Graph const& g )  
					{
						// reset color of all adjacent edges
						for ( auto ei = boost::out_edges( u, g ); ei.first !=  ei.second; ++ei.first ) {
							boost::put( boost::edge_color, g, *ei.first, 0 );
						}
					}
				// tree_edge
				template <typename Edge, typename Graph>	
					void tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 1 );
					}
				// non_tree_edge
				template <typename Edge, typename Graph>	
					void non_tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 2 );
					}
				// finish_vertex
				template <typename Vertex, typename Graph>	
					void finish_vertex( Vertex u, Graph const& g )  
					{
						auto it_end	= boost::get( vertex_cell, g, u )->incident_edge();
						auto it		= it_end->prev();
						while ( it != it_end ) {
							auto v = it->twin()->cell()->source_index();
							auto e = boost::edge( u, v, g ).first;
							if ( 2 == boost::get( boost::edge_color, g, e ) ) break;
							it = it->prev();
						}
						auto it2 = it;
						do {
							auto v = it2->twin()->cell()->source_index();
							auto e = boost::edge( u, v, g ).first;
							if ( 1 == boost::get( boost::edge_color, g, e ) ) {
								boost::add_edge( u, v, _SF );
							}
							it2 = it2->prev();
						} while ( it2 != it );
					}
		};

    }
	
    // ## struct        of planar_orderd_forest_solution
    
	template <
		typename CoordinateType     = double, 
		typename PointType			= boost::polygon::point_data<CoordinateType>,
		typename VoronoiDiagramType	= boost::polygon::voronoi_diagram<CoordinateType>	
	>
	class planar_orderd_forest_solution:
		// ## public inherit a boost undirected graph. It's associated properties is collected here.
		public 
			boost::adjacency_list
            <
				boost::vecS, 
				boost::vecS, 
				boost::undirectedS,

				boost::property<vertex_mother_t,		    std::size_t,		                // used to indicate the parent of the vertex 
				boost::property<vertex_cell_t,			    typename boost::polygon::voronoi_diagram<CoordinateType>::const_cell_iterator,
				boost::property<boost::vertex_color_t,      int,			
				boost::property<boost::vertex_distance_t,   CoordinateType			
				>>>>,

				boost::property<boost::edge_index_t,	int,				                // Used as a marker. 0: Triangulation Graph Edge, 1: MST Edge, 2: Boundary, 3: net
				boost::property<boost::edge_color_t,	int,				                // Used as a tempary marker for traversing
				boost::property<boost::edge_weight_t,	CoordinateType,		                // Used by kruskal's algorithm
				boost::property<boost::edge_flow_t,		int,			
				boost::property<boost::edge_capacity_t,	int,				                // Used as a recording of # of traces passing an edge
				boost::property<edge_point_t,	        std::array<std::list<int>, 2>		// Used as a recording of geometry passing points. The front side of the list is the lower id vertex
				>>>>>>
			> 
	{
		// ## sub-types
		public:
			
			using coordinate_type		= CoordinateType;
			using point_type			= PointType;
			using linestring_type		= boost::geometry::model::linestring<point_type>;
			using voronoi_diagram_type	= VoronoiDiagramType;
	
		// ## data members
        
		public:
			std::size_t			                    num_roots;
			std::size_t			                    num_leafs;
		private:
			std::vector<PointType>					_points;
			std::shared_ptr<VoronoiDiagramType>		_vd;

            std::vector< unsigned >                 _perturbs; 
			
		// ## function members 

        // ### constructor

		public:
			template <typename RootContainer, typename LeafContainer, typename NetContainer>
				planar_orderd_forest_solution( RootContainer&& roots, LeafContainer&& leafs, NetContainer&& nets );
		private:
			void construct_delaunay_trianguation( void );
            void mark_boundary_edges( void );
			template < typename NetContainer > 
                void mark_net_edges( NetContainer const& nets );

        // ### find_MST_by_comparable_length 
        
        public:
			void find_MST_by_comparable_length( void );
        private:
            template < typename MOMMap, typename IDXMap > 
                void mark_MST_edges( MOMMap v_mom_mp, IDXMap e_idx_mp );
	            
        // ### find_MST_by_heuristic_method 
        
        public:
			void find_MST_by_heuristic_method( void );
        private:
            template < typename IDXMap, typename WGTMap > 
                void set_edge_weight( IDXMap e_idx_mp, WGTMap e_wgt_mp );
            template < typename Vertex, typename IDXMap >
                double initial_weight( Vertex s, Vertex t, IDXMap e_idx_mp ); 

        // ### find_trees_by_random

        public:
	        void find_trees_by_random( void );

        // ### find_MaxST_by_comparable_length

        public:
	        void find_MaxST_by_comparable_length( void );

        // ### find_MinLevelST

        public:
	        void find_MinLevelST( void );
 
        // ### find_MST_by_A_star_heu_distance

        public:
	        void find_MST_by_A_star_heu_distance( void );

        // ### find_MST_by_target_dist_from_root

        public:
	        void find_MST_by_target_dist_from_root( void );

        // ### print

        public: 
			void    print( std::ostream& os = std::cout );

        // ### gnuplot

        public:
			template < typename GeometricPath = std::vector< linestring_type > >
				void gnuplot( GeometricPath const& gp = GeometricPath{} );
			template < typename GeometricPath = std::vector< linestring_type > >
				void gnuplot( GeometricPath const& gp, unsigned i );
			template < typename GeometricPath = std::vector< linestring_type > >
	            void gnuplot( GeometricPath const& gp, std::string GIF_name );
			template < typename GeometricPath = std::vector< linestring_type > >
	            void gnuplot_png( GeometricPath const& gp, std::string PNG_name, unsigned x, unsigned y );
			template < typename GeometricPath = std::vector< linestring_type > >
	            void gnuplot_svg( GeometricPath const& gp, std::string PNG_name, unsigned x, unsigned y );
			template < typename GeometricPath = std::vector< linestring_type > >
	            void gnuplot_svg( GeometricPath const& gp, std::string PNG_name, unsigned i, unsigned x, unsigned y );
	            
        private:
            auto write_point_data( std::ofstream& gnuplot_script_file );
            auto write_triangulation_data( std::ofstream& gnuplot_script_file );
            auto write_MST_data( std::ofstream& gnuplot_script_file );
            auto write_boundary_data( std::ofstream& gnuplot_script_file );
			template < typename GeometricPath >
	            auto write_path_data( std::ofstream& gnuplot_script_file, GeometricPath const& gp );
			template < typename GeometricPath >
	            auto write_path_data( std::ofstream& gnuplot_script_file, GeometricPath const& gp, unsigned i );

        // ### SA functions

        public:
			void    perturb( void );
			bool    switch_mother( std::size_t u );
			double  cost( void );
            void    write_perturb( std::ostream& os ) { for ( auto const& i : _perturbs ) { os << i << "\n"; } }

        // ### point
        
        public:
			auto    point( std::size_t u ) { return this->_points[ u ]; }

	};                                        

	// ## function def  of planar_orderd_forest_solution  
	
	// ### constructor 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template <typename RootContainer, typename LeafContainer, typename NetContainer>
	planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	planar_orderd_forest_solution( RootContainer&& roots, LeafContainer&& leafs, NetContainer&& nets ): 
		num_roots( roots.size() ), num_leafs( leafs.size() ) 
	{
		RootContainer	_roots( std::forward<RootContainer>( roots ) );
		LeafContainer	_leafs( std::forward<LeafContainer>( leafs ) );
		_points.insert( _points.end(), std::make_move_iterator( _roots.begin() ), std::make_move_iterator( _roots.end() ) );
		_points.insert( _points.end(), std::make_move_iterator( _leafs.begin() ), std::make_move_iterator( _leafs.end() ) );

		_vd = std::make_shared<VoronoiDiagramType>(); 
		boost::polygon::construct_voronoi( _points.begin(), _points.end(), &*_vd );
		construct_delaunay_trianguation();

	    mark_boundary_edges();
		mark_net_edges( nets );
	};

	// #### construct_delaunay_trianguation
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	construct_delaunay_trianguation( void )	
	{
		for ( auto cell_it = _vd->cells().begin(); cell_it != _vd->cells().end(); ++cell_it ) 
        {
			auto	u		    = cell_it->source_index();
			auto	edge_it	    = cell_it->incident_edge();

			cell_it->color( 1 );
			do 
            {
				if ( auto adj_cell_it = edge_it->twin()->cell(); adj_cell_it->color() == 0 ) 
                {
					add_edge( u, adj_cell_it->source_index(), *this );
				}
				edge_it = edge_it->prev();
			} while ( edge_it != cell_it->incident_edge() );

            put( vertex_cell, *this, u, cell_it );
		}
	}

	// #### mark_boundary_edges
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	mark_boundary_edges( void ) 
	{
		for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto v = source( *it, *this );
            auto u = target( *it, *this );
            if (    (v - u == 1 && v < num_roots) 
                 || (u - v == 1 && u < num_roots) )
            {
			    put( boost::edge_index, *this, *it, 2 );
            }
            /*
            if ( source( *it, *this ) < num_roots && 
                 target( *it, *this ) < num_roots     ) 
            {
			    put( boost::edge_index, *this, *it, 2 );
            }
            */
        }
    }

	// #### mark_net_edges

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename NetContainer >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	mark_net_edges( NetContainer const& nets )
	{
		for ( auto const& n: nets ) {
			auto [e, b] = add_edge( std::get<0>( n ), std::get<1>( n ), *this );
			put( boost::edge_index, *this, e, 3 );
		}
	}
	

	// ### find_MST_by_comparable_length
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MST_by_comparable_length( void ) 
	{
		auto v_mom_mp   = boost::get( vertex_mother             , *this );
		auto e_idx_mp   = boost::get( boost::edge_index         , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight        , *this ); 
        for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto e = *it;
            auto u = source( e, *this );
            auto v = target( e, *this );
            if      ( 0 == e_idx_mp[ e ] )
            {
				e_wgt_mp[ e ] = boost::geometry::comparable_distance( _points[ u ], _points[ v ] );
            }
            else if ( 2 == e_idx_mp[ e ] )
            {
                if ( 1 == u - v or 1 == v - u  ) 
                { 
                    e_wgt_mp[ e ] = 0; 
                }
                else
                {
                    e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
                }
            }
            else if ( 3 == e_idx_mp[ e ] )
            {
                e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
            }
        }
        boost::prim_minimum_spanning_tree( *this, v_mom_mp );
        mark_MST_edges( v_mom_mp, e_idx_mp );
    }

	// #### mark_MST_edges
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
    template < typename MOMMap, typename IDXMap >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	mark_MST_edges( MOMMap v_mom_mp, IDXMap e_idx_mp ) 
	{
        for ( auto [ vi, vi_end ] = vertices( *this ); vi != vi_end; ++vi )
        {
            if ( *vi < num_roots ) { v_mom_mp[ *vi ] = *vi; }
            if ( v_mom_mp[ *vi ] != *vi ) { e_idx_mp[ edge( *vi, v_mom_mp[ *vi ], *this ).first ] = 1; }
        }
    }


	// ### find_MST_by_heuristic_method
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MST_by_heuristic_method( void ) 
	{
		auto v_mom_mp   = boost::get( vertex_mother, *this );
		auto e_idx_mp   = boost::get( boost::edge_index, *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight, *this ); 

        set_edge_weight( e_idx_mp, e_wgt_mp );

        boost::prim_minimum_spanning_tree( *this, v_mom_mp );

        mark_MST_edges( v_mom_mp, e_idx_mp );
	}

	// #### set_edge_weight
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
    template < typename IDXMap, typename WGTMap >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	set_edge_weight( IDXMap e_idx_mp, WGTMap e_wgt_mp ) 
	{
        for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto e = *it;
            auto u = source( e, *this );
            auto v = target( e, *this );
            if      ( 0 == e_idx_mp[ e ] )
            {
                auto s = source( e, *this );
                auto t = target( e, *this );
                e_wgt_mp[ e ] = initial_weight( s, t, e_idx_mp );
            }
            else if ( 2 == e_idx_mp[ e ] )
            {
                if ( 1 == u - v or 1 == v - u  ) 
                { 
                    e_wgt_mp[ e ] = 0; 
                }
                else
                {
                    e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
                }
            }
            else if ( 3 == e_idx_mp[ e ] )
            {
                e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
            }
        }
    }

	// ##### initial_weight
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
    template < typename Vertex, typename IDXMap >
	double planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	initial_weight( Vertex s, Vertex t, IDXMap e_idx_mp ) 
	{
        auto N = boost::filtered_graph( *this, detail::is_net_edge( e_idx_mp ) );

        auto [s_it, s_it_end] = adjacent_vertices( s, N );
        auto [t_it, t_it_end] = adjacent_vertices( t, N );

        unsigned ds = s, dt = t;

        if ( s >= num_roots ) 
        {
            if ( s_it == s_it_end ) 
                return num_roots;
            else 
                ds = *s_it;
        }
        if ( t >= num_roots ) 
        {
            if ( t_it == t_it_end ) 
                return num_roots;
            else 
                dt = *t_it;
        }
        auto D          = (num_roots + ds - dt) % num_roots; 
        auto D_prime    = num_roots - D;
        return std::min( D, D_prime );
    }

	

	// ### find_trees_by_random
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_trees_by_random( void ) 
	{
		auto v_mom_mp   = boost::get( vertex_mother             , *this );
		auto e_idx_mp   = boost::get( boost::edge_index         , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight        , *this ); 

		auto engine = std::default_random_engine( 3241 );           // seed: 3241
		auto distr  = std::uniform_int_distribution<std::size_t>( 0, 10000 ); 
        for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto e = *it;
            auto u = source( e, *this );
            auto v = target( e, *this );
            if      ( 0 == e_idx_mp[ e ] )
            {
				e_wgt_mp[ e ] = distr( engine );
            }
            else if ( 2 == e_idx_mp[ e ] )
            {
                if ( 1 == u - v or 1 == v - u  ) 
                { 
                    e_wgt_mp[ e ] = 0; 
                }
                else
                {
                    e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
                }
            }
            else if ( 3 == e_idx_mp[ e ] )
            {
                e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
            }
        }

        boost::prim_minimum_spanning_tree( *this, v_mom_mp );
        mark_MST_edges( v_mom_mp, e_idx_mp );
	}

	// ### find_MaxST_by_comparable_length
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MaxST_by_comparable_length( void ) 
	{
		auto v_mom_mp   = boost::get( vertex_mother             , *this );
		auto e_idx_mp   = boost::get( boost::edge_index         , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight        , *this ); 
        for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto e = *it;
            auto u = source( e, *this );
            auto v = target( e, *this );
            if      ( 0 == e_idx_mp[ e ] )
            {
				e_wgt_mp[ e ] = 1 / boost::geometry::comparable_distance( _points[ u ], _points[ v ] );
            }
            else if ( 2 == e_idx_mp[ e ] )
            {
                if ( 1 == u - v or 1 == v - u  ) 
                { 
                    e_wgt_mp[ e ] = 0; 
                }
                else
                {
                    e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
                }
            }
            else if ( 3 == e_idx_mp[ e ] )
            {
                e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
            }
        }
        boost::prim_minimum_spanning_tree( *this, v_mom_mp );
        mark_MST_edges( v_mom_mp, e_idx_mp );
    }

	// ### find_MinLevelST
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MinLevelST( void ) 
	{
		auto v_mom_mp   = boost::get( vertex_mother         , *this );
		auto e_idx_mp   = boost::get( boost::edge_index     , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight    , *this ); 

        auto ColorMap   = std::vector<unsigned>( num_vertices(*this), 0 );       // 0: undiscovered, 1: discovered, 2: finished
        auto LevelMap   = std::vector<unsigned>( num_vertices(*this), 0 );
        auto Queue      = std::list<unsigned>();

        for ( unsigned i = 0; i != this->num_roots; ++i )
        {
            ColorMap[i] = 2;
            for ( auto [ei, ei_end] = out_edges( i, *this ); ei != ei_end; ++ei )
            {
                if ( 0 != e_idx_mp[ *ei ] ) continue;
                auto v = target( *ei, *this );
                if ( v >= this->num_roots )
                {
                    e_wgt_mp[ *ei ] = 1;
                    if ( ColorMap[v] == 0 )
                    {
                        LevelMap[v] = 1;
                        ColorMap[v] = 1;
                        Queue.push_back( v );
                    }
                }
            }
        }
        while ( !Queue.empty() )
        {
            auto u = Queue.front(); Queue.pop_front();
            for ( auto [ei, ei_end] = out_edges( u, *this ); ei != ei_end; ++ei )
            {
                if ( 0 != e_idx_mp[ *ei ] ) continue;
                auto v = target( *ei, *this );
                if ( ColorMap[v] == 2 ) continue;
                e_wgt_mp[ *ei ] = LevelMap[u] + 1;
                if ( ColorMap[v] == 0 )
                {
                    LevelMap[v] = LevelMap[u] + 1;
                    ColorMap[v] = 1;
                    Queue.push_back( v );
                }
            }
            ColorMap[u] = 2;
        }
        for ( auto [it, it_end] = edges( *this ); it != it_end; ++it ) 
        {
            auto e = *it;
            auto u = source( e, *this );
            auto v = target( e, *this );

            if ( 2 == e_idx_mp[ e ] )
            {
                if ( 1 == u - v or 1 == v - u  ) 
                { 
                    e_wgt_mp[ e ] = 0; 
                }
                else
                {
                    e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
                }
            }
            else if ( 3 == e_idx_mp[ e ] )
            {
                e_wgt_mp[ e ] = std::numeric_limits<CoordinateType>::max();
            }
        }
        boost::prim_minimum_spanning_tree( *this, v_mom_mp );
        mark_MST_edges( v_mom_mp, e_idx_mp );
    }

	// ### find_MST_by_A_star_heu_distance
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MST_by_A_star_heu_distance( void ) 
	{
		auto e_idx_mp   = boost::get( boost::edge_index     , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight    , *this ); 
		auto v_mom_mp   = boost::get( vertex_mother         , *this );

		auto F      = boost::filtered_graph( *this, detail::is_triangulation_edge(e_idx_mp) );
		auto nets   = boost::filtered_graph( *this, detail::is_net_edge(e_idx_mp) );

        for ( unsigned i = 0, i_end = num_vertices( *this ); i != i_end; ++i )
        {
            v_mom_mp[ i ] = i;
        }
        for ( auto [ei, ei_end] = edges( *this ); ei != ei_end; ++ei )
        {
		    e_wgt_mp[ *ei ] = boost::geometry::comparable_distance( _points[ source( *ei, *this ) ], _points[ target( *ei, *this ) ] );
        }

        using   DE          = std::tuple<unsigned, unsigned, double>;
        auto    cmp         = [&]( DE const& l, DE const& r ) { return std::get<2>(l) > std::get<2>(r); };
        auto    Q           = std::priority_queue< DE, std::vector<DE>, decltype(cmp) >( cmp );
        auto    ColorMap    = std::vector<unsigned>( num_vertices(*this), 0 );       // 0: undiscovered, 1: discovered, 2: finished
                                               
        for ( unsigned i = 0; i != this->num_roots; ++i )
        {
            for ( auto [ei, ei_end] = out_edges( i, F ); ei != ei_end; ++ei )
            {
                auto v = target( *ei, F );
                if ( v < this->num_roots )  { continue;         }
                if ( 0 == ColorMap[v] )     { ColorMap[v] = 1;  }
                Q.push( DE( i, v, e_wgt_mp[ *ei ] + e_wgt_mp[ *out_edges( v, nets ).first ] ) );
            }
            ColorMap[i] = 2;
        }
        while ( !Q.empty() )
        {
            auto e = Q.top(); Q.pop();
            if ( 1 == ColorMap[ std::get<1>(e) ] )
            {
                v_mom_mp[ std::get<1>(e) ] = std::get<0>(e);
            }
            for ( auto [ei, ei_end] = out_edges( std::get<1>(e), F ); ei != ei_end; ++ei )
            {
                auto v = target( *ei, F );
                if ( 2 == ColorMap[v] ) { continue;         }
                if ( 0 == ColorMap[v] ) { ColorMap[v] = 1;  }
                Q.push( DE( std::get<1>(e), v, e_wgt_mp[ *ei ] + e_wgt_mp[ *out_edges( v, nets ).first ] ) );
            }
            ColorMap[ std::get<1>(e) ] = 2;
        }

        mark_MST_edges( v_mom_mp, e_idx_mp );
    }

	// ### find_MST_by_target_dist_from_root
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	find_MST_by_target_dist_from_root( void ) 
	{
		auto e_idx_mp   = boost::get( boost::edge_index     , *this ); 
		auto e_wgt_mp   = boost::get( boost::edge_weight    , *this ); 
		auto v_mom_mp   = boost::get( vertex_mother         , *this );

		auto F      = boost::filtered_graph( *this, detail::is_triangulation_edge(e_idx_mp) );
		auto nets   = boost::filtered_graph( *this, detail::is_net_edge(e_idx_mp) );

        using   DE          = std::tuple<unsigned, unsigned, double>;
        auto    cmp         = [&]( DE const& l, DE const& r ) { return std::get<2>(l) > std::get<2>(r); };
        auto    Q           = std::priority_queue< DE, std::vector<DE>, decltype(cmp) >( cmp );
        auto    ColorMap    = std::vector<unsigned>( num_vertices(*this), 0 );       // 0: undiscovered, 1: discovered, 2: finished
        auto    RootMap     = std::vector<std::size_t>( num_vertices(*this) );       
                                               
        for ( unsigned i = 0, i_end = num_vertices( *this ); i != i_end; ++i )
        {
            v_mom_mp[ i ] = i;
            RootMap[ i ]  = i;
        }

        for ( unsigned i = 0; i != this->num_roots; ++i )
        {
            for ( auto [ei, ei_end] = out_edges( i, F ); ei != ei_end; ++ei )
            {
                std::size_t v = target( *ei, F );
                if ( v < this->num_roots )  { continue;         }
                if ( 0 == ColorMap[v] )     { ColorMap[v] = 1;  }
                Q.push( DE( i, v, initial_weight( RootMap[i], v, e_idx_mp ) ) );
            }
            ColorMap[i] = 2;
        }
        auto sz = Q.size();
        while ( !Q.empty() )
        {
            auto e = Q.top(); Q.pop();
            if ( 1 == ColorMap[ std::get<1>(e) ] )
            {
                v_mom_mp[ std::get<1>(e) ] = std::get<0>(e);
                RootMap[ std::get<1>(e) ] = std::get<0>(e);
            }
            for ( auto [ei, ei_end] = out_edges( std::get<1>(e), F ); ei != ei_end; ++ei )
            {
                std::size_t v = target( *ei, F );
                if ( 2 == ColorMap[v] ) { continue;        }
                if ( 0 == ColorMap[v] ) { ColorMap[v] = 1; }
                Q.push( DE( std::get<1>(e), v, initial_weight( RootMap[std::get<1>(e)], v, e_idx_mp ) ) );
            }
            ColorMap[ std::get<1>(e) ] = 2;
        }

        mark_MST_edges( v_mom_mp, e_idx_mp );
    }

	// ### print
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	print( std::ostream& os )
	{
		auto const& g = *this;
		unsigned col_width[5] = { 7, 35, 15, 13, 28 };

		auto header = 
			fmt::format( 
				"|-{0:-<{6}}-|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|\n" 
				"| {1: <{6}} | {2: ^{7}} | {3: <{8}} | {4: <{9}} | {5: <{10}} |\n"
				"|-{0:-<{6}}-|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|\n", 

				"", "vertex", "adjacency_list", "vertex_mother_t", "vertex_deg_t", "point get from vertex_cell_t",  
				col_width[0], col_width[1],	col_width[2], col_width[3], col_width[4]
			);

		os << header;

		for ( auto [it, it_end] = vertices( g ); it != it_end; ++it ) {

			auto [adj_it, adj_it_end] = adjacent_vertices( *it, g );
			std::vector<unsigned> adj_v( adj_it, adj_it_end );
			auto a1 = fmt::format( "{0}", adj_v );
			
			auto const& p = g._points[ get( vertex_cell, g, *it )->source_index() ];
			auto a2 = fmt::format( "({0: >10.2f}, {1: >10.2f})", p.x(), p.y() );

			auto row = 
				fmt::format(
					"| {0: <{5}} | {1: ^{6}} | {2: >{7}} | {3: >{8}} | {4: >{9}} |", 

					*it,			a1,				get( vertex_mother, g, *it ),	"x",	a2,  
					col_width[0],	col_width[1],	col_width[2],					col_width[3],				col_width[4]
				);

			os << row;	

			if (*it == (g.num_roots)) {
				os << fmt::format( "  <- num_roots: {}, num_leafs: {}", g.num_roots, g.num_leafs );
			}

			os << "\n";

		}


		unsigned col_width2[6] = { 10, 12, 12, 13, 11, 15 };

		auto header2 = 
			fmt::format( 
				"|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|-{0:-<{11}}-|-{0:-<{12}}-|\n" 
				"| {1: >{7}} | {2: >{8}} | {3: >{9}} | {4: >{10}} | {5: >{11}} | {6: >{12}} |\n" 
				"|-{0:-<{7}}-|-{0:-<{8}}-|-{0:-<{9}}-|-{0:-<{10}}-|-{0:-<{11}}-|-{0:-<{12}}-|\n" 

				, "", "edge", "edge_index_t", "edge_color_t", "edge_weight_t", "edge_flow_t", "edge_capacity_t"  
				, col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
			);

		os << header2;

		for (auto [it, it_end] = edges( g ); it != it_end; ++it) {

			auto edge = 
				fmt::format( "{{{0: >3}, {1: >3}}}", source( *it, g ), target( *it, g ) );

			auto w = get( boost::edge_weight, g, *it );
			auto weight = (w == std::numeric_limits<CoordinateType>::max())?  
				fmt::format( "{0: >{1}}", "max", col_width2[3] ): fmt::format( "{0: >{1}.2f}", w, col_width2[3] );
					

			auto row = 
				fmt::format( 
					"| {1: >{7}} | {2: >{8}} | {3: >{9}} | {4} | {5: >{11}} | {6: >{12}} |\n"

					, "", edge, get( boost::edge_index, g, *it ), get( boost::edge_color, g, *it ), weight, get( boost::edge_flow, g, *it ), get( boost::edge_capacity, g, *it ) 
					, col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
				);

			os << row;
		}	

		auto footage = 
			fmt::format( 
				"|={0:=<{1}}=|={0:=<{2}}=|={0:=<{3}}=|={0:=<{4}}=|={0:=<{5}}=|={0:=<{6}}=|\n" 

				, "", col_width2[0], col_width2[1], col_width2[2], col_width2[3], col_width2[4], col_width2[5]
			);

		os << footage;
	}

	// ### gnuplot 

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot( GeometricPath const& gp )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << fmt::format( ", {}    \\\n"  , path_plot_cmd             )
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
//      system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot( GeometricPath const& gp, unsigned i )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp, i );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << fmt::format( ", {}    \\\n"  , path_plot_cmd             )
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
        system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot( GeometricPath const& gp, std::string GIF_name )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << "set terminal gif animate delay 1 size 1000, 800 \n"
            << "set output '" << GIF_name << "'\n"
            << "set size ratio 1\n"
            << "unset xtics\n" 
            << "unset ytics\n" 
            << "n=" << gp.size() << "\n"
            << "do for [i=0:n] {\n"
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << ", '$paths' every :::::i using 1:2 with lines ls 4 notitle\n" 
            << "}\n"
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
        system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot_png( GeometricPath const& gp, std::string PNG_name, unsigned x, unsigned y )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << "set terminal png size " << x << "," << y << "\n"
            << "set output '" << PNG_name << ".png'\n"
            << "set size ratio 1\n"
            << "unset xtics\n" 
            << "unset ytics\n" 
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << fmt::format( ", {}    \\\n"  , path_plot_cmd             )
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
        system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot_svg( GeometricPath const& gp, std::string PNG_name, unsigned x, unsigned y )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << "set terminal svg size " << x << "," << y << "\n"
            << "set output '" << PNG_name << ".svg'\n"
            << "set size ratio 1\n"
            << "unset xtics\n" 
            << "unset ytics\n" 
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << fmt::format( ", {}    \\\n"  , path_plot_cmd             )
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
        system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	gnuplot_svg( GeometricPath const& gp, std::string PNG_name, unsigned i, unsigned x, unsigned y )
	{
        auto    gnuplot_script_file_name    = std::string( "sol.gp" );
        auto    gnuplot_script_file         = std::ofstream( gnuplot_script_file_name );

        auto    triangulation_plot_cmd      = write_triangulation_data( gnuplot_script_file );
        auto    boundary_plot_cmd           = write_boundary_data( gnuplot_script_file );
        auto    MST_plot_cmd                = write_MST_data( gnuplot_script_file );
        auto    path_plot_cmd               = write_path_data( gnuplot_script_file, gp, i );
        auto    points_plot_cmd             = write_point_data( gnuplot_script_file );

        gnuplot_script_file  
            << "set terminal svg size " << x << "," << y << "\n"
            << "set output '" << PNG_name << ".svg'\n"
            << "set size ratio 1\n"
            << "unset xtics\n" 
            << "unset ytics\n" 
            << fmt::format( "plot {} \\\n"  , triangulation_plot_cmd    )
            << fmt::format( ", {}    \\\n"  , boundary_plot_cmd         )
            << fmt::format( ", {}    \\\n"  , MST_plot_cmd              )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[0]        )
            << fmt::format( ", {}    \\\n"  , points_plot_cmd[1]        )
            << fmt::format( ", {}    \\\n"  , path_plot_cmd             )
        ;

        gnuplot_script_file.close();

        system( fmt::format( "gnuplot -p {}", gnuplot_script_file_name ).c_str() ); 
        system( fmt::format( "rm {}"        , gnuplot_script_file_name ).c_str() ); 
	}

    // #### write_point_data

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_point_data( std::ofstream& gnuplot_script_file )
    {
		gnuplot_script_file << "$points << EOD \n";

        for ( unsigned i = 0; i != _points.size(); ++i )
        {
			gnuplot_script_file << 
				fmt::format( 
					"{point_x} {point_y} {point_id} {pt_var} {ps_var} {lc_rgb_var} \n" 
					, "point_x"_a       = _points[i].x()
                    , "point_y"_a       = _points[i].y()
                    , "point_id"_a      = i 
					, "pt_var"_a        = 7
                    , "ps_var"_a        = 1
                    , "lc_rgb_var"_a    = "0x879cc4" 
				);
        }

		gnuplot_script_file << "EOD \n";
			
	    return std::array< std::string, 2 >
        {
			"'$points' using 1:2:4:5:6 with points pt variable ps variable lc rgb variable title '' ",
			"'$points' using 1:2:3 with labels offset char 0.5,0.5 title '' "
        };
    }

    // #### write_triangulation_data

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_triangulation_data( std::ofstream& gnuplot_script_file )
    {
		auto    e_idx_mp = 
                boost::get( boost::edge_index, *this ); 
		auto    triangulation =
                boost::filtered_graph( *this, detail::is_triangulation_edge( e_idx_mp ) );											

		gnuplot_script_file << "set style line 1 lt rgb 0xA9A9A9 \n";
        gnuplot_script_file << "$triangulation << EOD \n";
                
        for ( auto [ei, ei_end] = edges( triangulation ); ei != ei_end; ++ei )
        {
            auto const& ps = _points[ source( *ei, triangulation ) ]; 
            auto const& pt = _points[ target( *ei, triangulation ) ]; 
			gnuplot_script_file << 
                fmt::format( 
                    "{} {}\n" 
                    "{} {}\n" 
                    "\n"
                    , ps.x(), ps.y(), pt.x(), pt.y()   
                );
        }

        gnuplot_script_file << "EOD \n";
						
	    return std::string( "'$triangulation' using 1:2 with lines ls 1 title ''" );
    }

    // #### write_MST_data

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_MST_data( std::ofstream& gnuplot_script_file )
    {
		auto    e_idx_mp = 
                boost::get( boost::edge_index, *this ); 
		auto    MST =
                boost::filtered_graph( *this, detail::is_MST_edge( e_idx_mp ) );											

		gnuplot_script_file << "set style line 2 lt rgb 0xF26279 \n";
        gnuplot_script_file << "$MST << EOD \n";
                
        for ( auto [ei, ei_end] = edges( MST ); ei != ei_end; ++ei )
        {
            auto const& ps = _points[ source( *ei, MST ) ]; 
            auto const& pt = _points[ target( *ei, MST ) ]; 
			gnuplot_script_file << 
                fmt::format( 
                    "{} {}\n" 
                    "{} {}\n" 
                    "\n"
                    , ps.x(), ps.y(), pt.x(), pt.y()   
                );
        }

        gnuplot_script_file << "EOD \n";
						
	    return std::string( "'$MST' using 1:2 with lines ls 2 title ''" );
    }

    // #### write_boundary_data

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_boundary_data( std::ofstream& gnuplot_script_file )
    {
		auto    e_idx_mp = 
                boost::get( boost::edge_index, *this ); 
		auto    boundaries =
                boost::filtered_graph( *this, detail::is_boundary_edge( e_idx_mp ) );											

		gnuplot_script_file << "set style line 3 lt rgb 0x000000 \n";
        gnuplot_script_file << "$boundaries << EOD \n";
                
        for ( auto [ei, ei_end] = edges( boundaries ); ei != ei_end; ++ei )
        {
            auto const& ps = _points[ source( *ei, boundaries ) ]; 
            auto const& pt = _points[ target( *ei, boundaries ) ]; 
			gnuplot_script_file << 
                fmt::format( 
                    "{} {}\n" 
                    "{} {}\n" 
                    "\n"
                    , ps.x(), ps.y(), pt.x(), pt.y()   
                );
        }

        gnuplot_script_file << "EOD \n";
						
	    return std::string( "'$boundaries' using 1:2 with lines ls 3 title ''" );
    }

    // #### write_path_data

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_path_data( std::ofstream& gnuplot_script_file, GeometricPath const& gp )
    {
		gnuplot_script_file << "set style line 4 lt rgb 0x95D47A \n";
        gnuplot_script_file << "$paths << EOD \n";
                
        for ( auto it = std::begin( gp ), it_end = std::end( gp ); it != it_end; ++it )
        {
            for ( auto const& point : *it ) 
            {
                gnuplot_script_file << 
                    fmt::format( 
                        "{} {}\n" 
                        , boost::geometry::get<0>(point)
                        , boost::geometry::get<1>(point)
                    );
            }
            gnuplot_script_file << "\n";
        }

        gnuplot_script_file << "EOD \n";
						
	    return std::string( "'$paths' using 1:2 with lines ls 4 notitle" );
    }

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	template < typename GeometricPath >
	auto planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	write_path_data( std::ofstream& gnuplot_script_file, GeometricPath const& gp, unsigned i )
    {
		gnuplot_script_file << "set style line 4 lt rgb 0x95D47A \n";
        gnuplot_script_file << "$paths << EOD \n";
                
        for ( auto it = std::begin( gp ), it_end = std::end( gp ); it != it_end; ++it )
        {
            for ( auto const& point : *it ) 
            {
                gnuplot_script_file << 
                    fmt::format( 
                        "{} {}\n" 
                        , boost::geometry::get<0>(point)
                        , boost::geometry::get<1>(point)
                    );
            }
            gnuplot_script_file << "\n";
        }

        gnuplot_script_file << "EOD \n";
						
	    return fmt::format( "'$paths' every :::::{} using 1:2 with lines ls 4 notitle", i );
    }


	// ### perturb 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	void planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	perturb( void )
	{
		std::default_random_engine 
			engine( std::chrono::system_clock::now().time_since_epoch().count() );
		std::uniform_int_distribution<std::size_t> 
			d( num_roots, num_leafs + num_roots - 1 );

        while (true)  
        {
            auto n = d( engine );
            _perturbs.push_back( n );
            if ( switch_mother( n ) ) { break; }
        }  
		//while ( !switch_mother( d( engine ) ) );
	}

	// ### switch_mother

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	bool planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	switch_mother( std::size_t u )
	{
		boost::filtered_graph											
			tri_edges( *this, detail::is_triangulation_edge( get( boost::edge_index, *this ) ) );

		auto const mu = get( vertex_mother, *this, u );
		auto edge_it_end = get( vertex_cell, *this, u )->incident_edge();
		while ( edge_it_end->twin()->cell()->source_index() != mu ) { edge_it_end = edge_it_end->prev(); }
		
		for ( auto edge_it = edge_it_end->prev(); edge_it != edge_it_end; edge_it = edge_it->prev() ) 
		{
			auto const x = edge_it->twin()->cell()->source_index();
			auto const [e1, b1] = edge( u, x, tri_edges );
			if ( not b1 ) continue;
			
			for ( auto v = x, w = get( vertex_mother, *this, v ); w != u; v = w, w = get( vertex_mother, *this, v ) )
			{
				if ( w == v )	// there is no loop
				{
					auto [e2, b2] = edge( u, mu, *this );
					put( boost::edge_index, *this, e2, 0 );
					put( boost::edge_index, *this, e1, 1 );
					put( vertex_mother, *this, u, x );
					return true;	
				}
			}
		} 
		return false;
	}

    // ### cost 
	
	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	double planar_orderd_forest_solution< CoordinateType, PointType, VoronoiDiagramType >::
	cost( void )
	{
        unsigned cross_number = 0;
		boost::filtered_graph F( *this, detail::is_MST_edge( get( boost::edge_index, *this ) ) );
		boost::filtered_graph nets( *this, detail::is_net_edge( get( boost::edge_index, *this ) ) );
        
        /* version 1

		boost::depth_first_search( F, visitor( detail::accumulate_cross_number( num_roots, num_vertices( F ), nets, cross_number ) ) );

        */
        /* version 2

        // BFS traversal to obtain ordered foreset
        boost::adjacency_list<> 
            ordered_forest;
		for ( std::size_t i = 0; i != num_roots; ++i ) {
			boost::breadth_first_search( 
				F, i , 
	            visitor( detail::build_sketchable_forest( ordered_forest ) )
			);
		}
        // DFS traversal to obtain dest_order and compute cross_number
		boost::depth_first_search( 
			ordered_forest, 
			visitor( detail::compute_dest_ordering_and_cross_number( nets, num_vertices( ordered_forest ), num_roots, cross_number ) )
		);
        */

        // version 3

        // BFS traversal to obtain ordered foreset
        boost::adjacency_list<> 
            ordered_forest;
		for ( std::size_t i = 0; i != num_roots; ++i ) {
			boost::breadth_first_search( 
				F, i , 
	            visitor( detail::build_sketchable_forest( ordered_forest ) )
			);
		}
        // DFS traversal to obtain dest_order and compute cross_number
		boost::depth_first_search( 
			ordered_forest, 
			visitor( detail::compute_dest_ordering_and_cross_number2( nets, num_vertices( ordered_forest ), num_roots, cross_number ) )
		);
        //

		return cross_number;
	}



	// # circular_frame
	
    namespace detail
    {
		// ## build_directed_forest
		
		template < typename DiGraph >
		class build_directed_forest: public boost::default_bfs_visitor 
		{
			// constructor 
			public:
				build_directed_forest( DiGraph& sf, unsigned n ): _dg( sf ), _num_roots( n ) {}
			// data member
			private:
				DiGraph&    _dg;
                unsigned    _num_roots;
			// events
			public:	
				// examine_vertex
				template <typename Vertex, typename Graph>	
					void examine_vertex( Vertex u, Graph const& g )  
					{
						// reset color of all adjacent edges
						for ( auto ei = boost::out_edges( u, g ); ei.first !=  ei.second; ++ei.first ) {
							boost::put( boost::edge_color, g, *ei.first, 0 );
						}
					}
				// tree_edge
				template <typename Edge, typename Graph>	
					void tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 1 );
					}
				// non_tree_edge
				template <typename Edge, typename Graph>	
					void non_tree_edge( Edge e, Graph const& g )  
					{
						boost::put( boost::edge_color, g, e, 2 );
					}
				// finish_vertex
				template <typename Vertex, typename Graph>	
					void finish_vertex( Vertex u, Graph const& g )  
					{
						auto it_end	= boost::get( vertex_cell, g, u )->incident_edge();
						auto it		= it_end->prev();
                        if ( u < _num_roots )
                        {
                            auto v = (_num_roots + u - 1) % _num_roots;
                            while ( v != it->twin()->cell()->source_index() ) { it = it->prev(); }
                        }
                        else
                        {
                            while ( it != it_end ) {
                                auto v = it->twin()->cell()->source_index();
                                auto e = boost::edge( u, v, g ).first;
                                if ( 2 == boost::get( boost::edge_color, g, e ) ) break;
                                it = it->prev();
                            }
                        }
						auto it2 = it;
						do {
							auto v = it2->twin()->cell()->source_index();
							auto e = boost::edge( u, v, g ).first;
							if ( 1 == boost::get( boost::edge_color, g, e ) ) {
								boost::add_edge( u, v, _dg );
							}
							it2 = it2->prev();
						} while ( it2 != it );
					}
		};

        // ## determine_net_priority_and_copmute_dest_ordering
        
		template < typename CircularFrame, typename Net, typename DestOrder >
		struct determine_net_priority_and_copmute_dest_ordering: public boost::default_dfs_visitor 
		{
            CircularFrame&  C;
            Net             N;
            int             B;      // number of roots on the boudary
            int             R;      // current root
            DestOrder&      D;

            determine_net_priority_and_copmute_dest_ordering( CircularFrame& c, Net& n, unsigned b, DestOrder& d ):
                C( c ), 
                N( n ), 
                B( b ),
                D( d )
            {}

			// start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
			// discover_vertex
			template < typename Vertex, typename Graph >
			void discover_vertex( Vertex u, Graph const& g ) 
			{
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    add_edge( u, *beg, C.nets ); 
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    D[ u ][ 0 ] = ( Dcw < Dccw )? Dcw : -Dccw;
                    D[ u ][ 1 ] = R;
                }
                else
                {
                    D[ u ][ 0 ] = 0;
                    D[ u ][ 1 ] = R;
                } 
            }
        };
        
        // ## compute_net_priority_and_dest_ordering
        
		template < typename CircularFrame, typename Net, typename DestOrder >
		struct compute_net_priority_and_dest_ordering: public boost::default_dfs_visitor 
		{
            using T = std::vector<std::vector<std::size_t>>;
            CircularFrame&  C;
            Net             N;
            int             B;      // number of roots on the boudary
            int             R;      // current root
            DestOrder&      D;

            T               net_order;
            unsigned        level;

            compute_net_priority_and_dest_ordering( CircularFrame& c, Net& n, unsigned b, DestOrder& d ):
                C( c ), 
                N( n ), 
                B( b ),
                D( d ),
                level( 0 )
            {}

			// start_vertex
			template < typename Vertex, typename Graph >
			void start_vertex( Vertex u, Graph const& g ) 
			{
				R = u;
			}
			// discover_vertex
			template < typename Vertex, typename Graph >
			void discover_vertex( Vertex u, Graph const& g ) 
			{
                compute_dest_ordering( u );
                if ( level < net_order.size() )
                {
                    net_order[ level ].push_back( u );
                }
                else 
                {
                    net_order.emplace_back( 1, u );         // push back a vector containing a u
                }
                ++level;
            }
			// finish_vertex
			template < typename Vertex, typename Graph >
			void finish_vertex( Vertex u, Graph const& g ) 
			{
                --level;
                if ( u == B - 1 )   // it's the last root
                {
                    for ( unsigned i = 1; i != net_order.size(); ++i )
                    {
                        for ( auto const& v: net_order[i] )
                        {
                            auto [ beg, end ] = adjacent_vertices( v, N );
                            add_edge( v, *beg, C.nets ); 
                        }
                    }
                }
			}
            // compput_dest_ordering
			template < typename Vertex >
            void compute_dest_ordering( Vertex u )
            {
                auto [ beg, end ] = adjacent_vertices( u, N );
                if ( u >= B && beg != end )  
                {
                    int Dcw  = ( B + R - *beg ) %  B;
                    int Dccw = B - Dcw;
                    D[ u ][ 0 ] = ( Dcw < Dccw )? Dcw : -Dccw;
                    D[ u ][ 1 ] = R;
                }
                else
                {
                    D[ u ][ 0 ] = 0;
                    D[ u ][ 1 ] = R;
                } 
            }
        };
        
    }
	template < typename CircularFrameSolution >
	class circular_frame :
		public syc::topology::model::v_01::sketchable_forest< 
			typename CircularFrameSolution::coordinate_type, 
			typename CircularFrameSolution::point_type 
		>
	{
		// ## declaration

		private:
			using coordinate_type	= typename CircularFrameSolution::coordinate_type;
			using point_type		= typename CircularFrameSolution::point_type;
			using linestring_type	= typename CircularFrameSolution::linestring_type;
			using base				= syc::topology::model::v_01::sketchable_forest< point_type, coordinate_type >;

		// ## member data
		
		private:
			CircularFrameSolution&		_sol;

		// ## member functions 

		public:
			circular_frame( CircularFrameSolution& sol );
			auto the_basic_routing_algorithm( void );
			auto t_escape( void );              // origin                                                   
			auto t_escape2( void );             // 修掉 root 不符合 t_escape 概念的地方                     
			auto t_escape3( void );             // 修掉 target 落在不同 slice 時不符合 t_escape 概念的地方  
			auto t_escape4( void );             // 改成 level first                                        
			auto realization( double spacing = std::numeric_limits<double>::max() );
            auto get_nets( void );
		private:
			bool clockwise_depth_first_route( typename base::topo_descriptor& h, std::size_t const u );
            template < typename T1, typename T2 >
                auto anticlockwise_route_through( T1 h, T2 T );
            template < typename T1, typename T2 >
	            auto clockwise_route_through( T1 h, T2 T );
	        auto find_entrances_to_diff_slice( typename base::topo_descriptor h );
	        auto anticlockwise_route_through( typename base::topo_descriptor h, typename base::topo_descriptor t );
	        auto clockwise_route_through( typename base::topo_descriptor h, typename base::topo_descriptor t );
            
	};

	// ## Methods of circular_frame
	
	// ### constructor
	
	template < typename CircularFrameSolution >
	circular_frame< CircularFrameSolution >::
	circular_frame( CircularFrameSolution& sol ): _sol( sol )
	{
		boost::filtered_graph											
			spannig_tree( _sol, detail::is_MST_edge( get( boost::edge_index, _sol ) ) );

		for ( std::size_t i = 0; i != _sol.num_roots; ++i ) {
			boost::breadth_first_search( 
				spannig_tree, i , 
				visitor( detail::build_directed_forest( *this, _sol.num_roots ) ) 
			);
		}
		this->make_sketch();
	}

	// ### the_basic_routing_algorithm
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	the_basic_routing_algorithm( void )
	{
		unsigned failure_cnt = 0; 

		boost::filtered_graph											
			nets( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );

        unsigned loop = 0;
		for ( auto [ ei, ei_end ] = edges( nets ); ei != ei_end; ++ei ) 
		{
			auto source_vertex = source( *ei, nets );					
			auto target_vertex = target( *ei, nets );					
			add_edge( source_vertex, target_vertex, this->nets );		// very ugly
            
//         std::cout << loop++ << std::endl;
//         if ( loop > 82 ) continue; 

			auto head = this->topo_instances( source_vertex )[ 0 ];	
			this->reset_slice_indexes();		// All the slices are not yet visited.
			if ( !clockwise_depth_first_route( head, target_vertex ) ) { ++failure_cnt; } 
		}	

		return failure_cnt;
	}

	// ### t_escape

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape( void )
	{
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::determine_net_priority_and_copmute_dest_ordering( *this, net, _sol.num_roots, dest_order ) ) 
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto T    = this->topo_instances( net_target );	
            auto input_it = std::find_if( std::begin( T ), std::end( T ), 
                [&]( auto t )
                {
                   return this->slice( h ) != this->slice( t ); 
                } 
            );

            this->reset_slice_indexes();            // All the slices are not yet visited.

            if ( input_it != std::end( T ) )
            {
                if ( !clockwise_depth_first_route( h, net_target ) ) { ++num_failed_nets; } 
            }
            else if ( dest_order[ net_source ][ 0 ] < 0 )       // should escape to the left of the current roots
            {
                auto anticlockwise_route_through = [&] ( auto& h, auto const& T )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; ++p != ori_mark; ) {
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                    for ( auto p = ori_mark; ++p != h; ) { 
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 0 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                };
                while ( !anticlockwise_route_through( h, T ) );
            }
            else
            {
                auto clockwise_route_through = [&] ( auto& h, auto const& T )
                {
                    auto m = this->slice( h );
                    auto ori_mark = m->end();
                    for ( auto p = h; --p != ori_mark; ) {
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );															
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                    for ( auto p = ori_mark; --p != h; ) { 
                        for ( auto const t: T )
                        {
                            if ( p == t )
                            {
                                h = this->make_slice( h, p );
                                return true;
                            } 
                        }
                        if ( this->attribution( p ) == 1 ) { 
                            auto e = this->ordering_pair( p );												
                            if ( this->slice( e ) == m ) { 													
                                h = this->make_slice( h, p );																	
                                return false;
                            }
                        }
                    }
                };
                while ( !clockwise_route_through( h, T ) );
            }
		}

		return num_failed_nets;
	}

	// ### t_escape2

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape2( void )
	{
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::determine_net_priority_and_copmute_dest_ordering( *this, net, _sol.num_roots, dest_order ) )   
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto T = this->topo_instances( net_target );	
            auto input_it = std::find_if( std::rbegin( T ), std::rend( T ), 
                [&]( auto t )
                {
                   return this->slice( h ) == this->slice( t ); 
                } 
            );

            this->reset_slice_indexes();            // all the slices are not yet visited.

            if ( input_it == std::rend( T ) )        // all the topoinstance of the net target are in deferent slices
            {
                if ( !clockwise_depth_first_route( h, net_target ) ) { ++num_failed_nets; } 
            }
            else if ( dest_order[ net_source ][ 0 ] < 0 )       // should escape to the left of the current roots
            {
                do
                {
                    h = anticlockwise_route_through( h, T );
                } while ( this->attribution( h ) != 2 );
            }
            else if ( dest_order[ net_source ][ 0 ] > 0 )       // should escape to the right of the current roots
            {
                do
                {
                    h = clockwise_route_through( h, T );
                } while ( this->attribution( h ) != 2 );
            }
            else
            {
                if ( input_it == std::rbegin( T ) )
                {
                    do
                    {
                        h = anticlockwise_route_through( h, T );
                    } while ( this->attribution( h ) != 2 );
                }
                else
                {
                    do
                    {
                        h = clockwise_route_through( h, T );
                    } while ( this->attribution( h ) != 2 );
                }
            }
		}

		return num_failed_nets;
    }

	// ### t_escape3

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape3( void )
	{
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::determine_net_priority_and_copmute_dest_ordering( *this, net, _sol.num_roots, dest_order ) )   
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto T = this->topo_instances( net_target );	

            this->reset_slice_indexes();                                        // all the slices are not yet visited.

            auto it = std::find_if( std::begin( T ), std::end( T ),             // find the first topo-instance in T that is in the same slice
                [&]( auto t )
                {
                   return this->slice( h ) == this->slice( t ); 
                } 
            );

            if ( it != std::end( T ) )                                          // if one of the topo-instances in T is in the same slice
            {
                if ( this->shortest_direction( h, *it ) )                         // if the shorest distance to t is going in a CCW direction
                {
                    h = anticlockwise_route_through( h, *it );
                }
                else                                                            // else if the shorest distance to t is going in a CW direction
                {
                    h = clockwise_route_through( h, *it );
                }
            }
            else
            {
                ++num_failed_nets;                                              // assume it's a failed net
                auto E = find_entrances_to_diff_slice( h );
                for ( auto e : E )
                {
                    if ( this->shortest_direction( h, e ) )                     // if the shorest distance to t is going in a CCW direction
                    {
                        h = anticlockwise_route_through( h, e );
                    }
                    else                                                        // else if the shorest distance to t is going in a CW direction
                    {
                        h = clockwise_route_through( h, e );
                    }

                    if ( clockwise_depth_first_route( h, net_target ) ) 
                    { 
                        --num_failed_nets;                                      // the net is sucessfully routed
                        break; 
                    } 
                    else
                    {
                        while ( this->attribution( h ) != 2 )
                        {
                            h = this->free_slice( h );
                        }
                    }
                }
            }
		}

		return num_failed_nets;
    }

	// ### t_escape4

	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	t_escape4( void )
	{
		unsigned num_failed_nets = 0;

        auto net        = boost::filtered_graph( _sol, detail::is_net_edge( get( boost::edge_index, _sol ) ) );
        auto dest_order = std::vector< std::array< int, 2 > >( num_vertices( *this ) ); // [ offset, currenet_root ]
		boost::depth_first_search( *this, 
            visitor( detail::compute_net_priority_and_dest_ordering( *this, net, _sol.num_roots, dest_order ) )   
        );

		for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei ) 
		{
			auto net_source = source( *ei, this->nets );
			auto net_target = target( *ei, this->nets );

            auto h = this->topo_instances( net_source )[ 0 ];	
            auto T = this->topo_instances( net_target );	

            this->reset_slice_indexes();                                        // all the slices are not yet visited.

            auto it = std::find_if( std::begin( T ), std::end( T ),             // find the first topo-instance in T that is in the same slice
                [&]( auto t )
                {
                   return this->slice( h ) == this->slice( t ); 
                } 
            );

            if ( it != std::end( T ) )                                          // if one of the topo-instances in T is in the same slice
            {
                if ( this->shortest_direction( h, *it ) )                         // if the shorest distance to t is going in a CCW direction
                {
                    h = anticlockwise_route_through( h, *it );
                }
                else                                                            // else if the shorest distance to t is going in a CW direction
                {
                    h = clockwise_route_through( h, *it );
                }
            }
            else
            {
                ++num_failed_nets;                                              // assume it's a failed net
                auto E = find_entrances_to_diff_slice( h );
                for ( auto e : E )
                {
                    if ( this->shortest_direction( h, e ) )                     // if the shorest distance to t is going in a CCW direction
                    {
                        h = anticlockwise_route_through( h, e );
                    }
                    else                                                        // else if the shorest distance to t is going in a CW direction
                    {
                        h = clockwise_route_through( h, e );
                    }

                    if ( clockwise_depth_first_route( h, net_target ) ) 
                    { 
                        --num_failed_nets;                                      // the net is sucessfully routed
                        break; 
                    } 
                    else
                    {
                        while ( this->attribution( h ) != 2 )
                        {
                            h = this->free_slice( h );
                        }
                    }
                }
            }
		}

		return num_failed_nets;
    }

	// ### realization
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	realization( double spacing ) 
	{
		auto R      = this->results();															    // every net should have a head before calling this method
		auto TCPs   = std::vector<std::vector<std::tuple<unsigned, unsigned, int, int>>>();	        // triangulation crossing paths
		for ( auto& path : R )  
		{						
			// step1: relax each of the paths
            auto dirty_patch = path.back();
			for ( bool relax = true; relax == true; ) {
				relax = false;
				auto it1 = path.begin();
				while ( true ) {
					auto it2 = std::next( it1 );
					if ( it2 == path.end() )	break;
					if ( it2->second == 0 )		break;
					if ( it2->second > 0 ) {
						auto cell_it = boost::get( vertex_cell, _sol, it1->first );
						auto vor_edge_it = cell_it->incident_edge();
						while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->next();	
						for ( ;    std::next( it2 )->first == vor_edge_it->next()->twin()->cell()->source_index()
								&& std::next( it2 )->second >= 0 && it1->second >=0 ;
							  vor_edge_it = vor_edge_it->next() ) {
							it2 = path.erase( it2 );
							relax = true;
						}
					}
					else if ( it2->second < 0 ) {
						auto cell_it = boost::get( vertex_cell, _sol, it1->first );
						auto vor_edge_it = cell_it->incident_edge();
						while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) vor_edge_it = vor_edge_it->prev();
						for ( ;    std::next( it2 )->first == vor_edge_it->prev()->twin()->cell()->source_index() 
								&& std::next( it2 )->second <= 0 && it1->second <=0 ;
							  vor_edge_it = vor_edge_it->prev() ) {
							it2 = path.erase( it2 );
							relax = true;
						}
					}
					it1 = it2;
				}
			}
            if ( path.back() != dirty_patch ) path.push_back( dirty_patch );
			// step2: realize each of the paths to a triangulation crossing paths
			std::vector<std::tuple<unsigned, unsigned, int, int>> X;
			if ( path.empty() )
			{
				TCPs.push_back( std::move( X ) );
			} 
			else
			{
				X.emplace_back( path.front().first, path.front().first, 0, 0 );
				auto	it1			= path.begin(), 
						it2			= std::next( it1 ),
						it3			= std::next( it2 );
				auto	vor_edge_it	= boost::get( vertex_cell, _sol, it1->first )->incident_edge();
				while ( vor_edge_it->twin()->cell()->source_index() != it2->first ) 
				{
					vor_edge_it = vor_edge_it->next();
				}
				for	( ;it2->second != 0; ++it1, ++it2, ++it3 )
				{
					vor_edge_it = vor_edge_it->twin();
					if ( it2->second > 0 ) 
					{ 
						if ( it1->second < 0 )
						{
							X.emplace_back( it2->first, it1->first, it2->second, it1->second );
						}
						vor_edge_it = vor_edge_it->prev();
						for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
							  v != it3->first; 
							  vor_edge_it = vor_edge_it->prev(), v = vor_edge_it->twin()->cell()->source_index() )  
						{
							X.emplace_back( it2->first, v, it2->second, 0 );
						}
					} 
					else if ( it2->second < 0 ) 
					{
						if ( it1->second > 0 )
						{
							X.emplace_back( it2->first, it1->first, it2->second, it1->second );
						}
						vor_edge_it = vor_edge_it->next();
						for ( auto v = vor_edge_it->twin()->cell()->source_index(); 
							  v != it3->first; 
							  vor_edge_it = vor_edge_it->next(), v = vor_edge_it->twin()->cell()->source_index() )  
						{
							X.emplace_back( it2->first, v, it2->second, 0 );
						}
					}
				}
				X.emplace_back( path.back().first, path.back().first, 0, 0 );
				TCPs.push_back( std::move( X ) );
			}
		}
        // step3: sort on edge
        auto GPds   = std::vector<std::vector<std::pair<typename std::list<int>::iterator,int>>>();
		for ( auto const& X : TCPs ) 
		{		
            auto gpd = std::vector<std::pair<typename std::list<int>::iterator,int>>();
			for ( auto const& t : X ) 
            { 
                auto [u, v, w, w2]  = t;
                auto n              = std::abs( w  );
                auto n2             = std::abs( w2 );
                if ( u == v ) { ; }
                else if ( u < v )    
                {
                    auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                    auto it = E[0].end(), it_end = it;
                    while ( ++it != it_end && *it < n ) { ; }
                    auto pos = E[0].insert( it, n );
                    gpd.emplace_back( pos, 0 );
                }
                else            
                {
                    if ( w2 == 0 )
                    {
                        auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                        auto it = E[1].end(), it_end = it;
                        while ( ++it != it_end && *it < n ) { ; }
                        auto pos = E[1].insert( it, n );
                        gpd.emplace_back( pos, 1 );
                    }
                    else
                    {
                        auto& E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
                        auto it = E[0].end(), it_end = it;
                        while ( ++it != it_end && *it < n2 ) { ; }
                        auto pos = E[0].insert( it, n2 );
                        gpd.emplace_back( pos, 0 );
                    }
                }
            }
            GPds.push_back( std::move( gpd ) );
        }
        // step4: transform to geometrical paths
		auto GPs    = std::vector<linestring_type>();			            // geometrical paths
        for ( unsigned i = 0; i != TCPs.size(); ++i )
        {
			linestring_type gpath;
            for ( unsigned j = 0; j != TCPs[i].size(); ++j )
            {
                auto const& [u, v, w, w2] = TCPs[i][j];
				auto const& p0 = ( w2 != 0 && u > v )? _sol.point( v ): _sol.point( u );
				auto const& p1 = ( w2 != 0 && u > v )? _sol.point( u ): _sol.point( v );
				auto&       E = boost::get( edge_point, _sol, boost::edge( u, v, _sol ).first );
				if ( u == v )
				{
				    auto const& p0 = _sol.point( u );
					auto const& x = boost::geometry::get<0>(p0);
					auto const& y = boost::geometry::get<1>(p0);
					boost::geometry::append( gpath, point_type( x, y ) );
				}
				else 
				{
                    auto d       = (double)std::distance( E[GPds[i][j-1].second].begin(), GPds[i][j-1].first );
                    auto frac    = std::min( (double)1 / (E[0].size() + E[1].size() + 1), spacing / boost::geometry::distance( p0, p1 ) );
                    auto offsetx = ( boost::geometry::get<0>(p0) - boost::geometry::get<0>(p1) ) * frac * (d + 1);
                    auto offsety = ( boost::geometry::get<1>(p0) - boost::geometry::get<1>(p1) ) * frac * (d + 1);
    				auto x       = boost::geometry::get<0>(p0) - offsetx;
    				auto y       = boost::geometry::get<1>(p0) - offsety;

					boost::geometry::append( gpath, point_type( x, y ) );
				}
			}
			GPs.push_back( std::move( gpath ) );
        }
		return GPs;
	}

	// ### get_nets
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	get_nets( void ) 
	{
        auto r = std::vector< std::array< std::size_t, 2 > >();
        for ( auto [ ei, ei_end ] = edges( this->nets ); ei != ei_end; ++ei )
        {
            auto u = source( *ei, this->nets );
            auto v = target( *ei, this->nets );
            r.emplace_back( std::array< std::size_t, 2 >{ std::min( u , v ), std::max( u, v ) } );
        } 
        return r;
    }

	// ### clockwise_depth_first_route			
	
	template < typename CircularFrameSolution >
	bool circular_frame< CircularFrameSolution >::
	clockwise_depth_first_route( typename base::topo_descriptor& h, std::size_t const u ) 
	{
		// Clockwise seach for the first target vertex t of h that is within the same slice of h.
		// If sucessfully finds a t, route from h to t and move on to the next source vertex.
		auto m = this->slice( h ); 																		// Let m = the pointer to slice where h is in.
		m->index = 1; 																				// Set the color of slice m gray
		for ( auto const& t : this->topo_instances( u ) ) { 
			auto n = this->slice( t ); 																	// Let n = the pointer to slice where target pointed by t is in.
			if ( m == n ) {
				this->make_slice( h, t );	
				return true; 																		// Successfully find a topological routing sketch
			}
		}
		// EVENT POINT: cannot route from h to any of its targets in the current slice
		
		// Anticlockwise search for any candidate edges to go through
		std::stack<typename base::topo_descriptor> E;
		auto ori_mark = this->slice( h )->end();
		for ( auto p = h; ++p != ori_mark; ) {
			if ( this->attribution( p ) != 2 ) { 
				auto e = this->ordering_pair( p );															// Let e be the related edge of the edge pointed by p
				if ( 0 == this->slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
					E.push( p );																		// the edge pointed by pi is chosen as a candidate.
				}
			}
		}
		for ( auto p = ori_mark; ++p != h; ) { 
			if ( this->attribution( p ) != 2 ) { 
				auto e = this->ordering_pair( p );															// Let e be the related edge of the edge pointed by p
				if ( 0 == this->slice( e )->index ) { 														// If the color of the slice is white, i.e. the slice has not been visited,
					E.push( p );																		// the edge pointed by pi is chosen as a candidate.
				}
			}
		}
		while ( !E.empty() ) {
			auto inArrow = E.top();
			h = this->make_slice( h, inArrow );					
			if ( clockwise_depth_first_route( h, u ) ) { 
				return true; 
			}
			E.pop();
		}

		// EVENT POINT: In the current slice, there is no target to route to, and no edge go through. I.e Have to turn around, back to the previous slice.
		
		h = this->free_slice( h );

		return false;
	}

	// ### anticlockwise_route_through			
	
	template < typename CircularFrameSolution >
    template < typename T1, typename T2 >
	auto circular_frame< CircularFrameSolution >::
	anticlockwise_route_through( T1 h, T2 T ) 
	{
        auto m = this->slice( h );
        auto ori_mark = m->end();
        for ( auto p = h; ++p != ori_mark; ) {
            for ( auto const t: T )
            {
                if ( p == t )
                {
                    return this->make_slice( h, p );
                } 
            }
            if ( this->attribution( p ) == 0 ) 
            { 
                auto e = this->ordering_pair( p );															
                if ( this->slice( e ) == m ) 
                { 													
                    return this->make_slice( h, p );																	
                }
            }
        }
        for ( auto p = ori_mark; ++p != h; ) { 
            for ( auto const t: T )
            {
                if ( p == t )
                {
                    return this->make_slice( h, p );
                } 
            }
            if ( this->attribution( p ) == 0 ) { 
                auto e = this->ordering_pair( p );												
                if ( this->slice( e ) == m ) 
                { 													
                    return this->make_slice( h, p );																	
                }
            }
        }
    }
    
	// ### clockwise_route_through			
	
	template < typename CircularFrameSolution >
    template < typename T1, typename T2 >
	auto circular_frame< CircularFrameSolution >::
	clockwise_route_through( T1 h, T2 T ) 
	{
        auto m = this->slice( h );
        auto ori_mark = m->end();
        for ( auto p = h; --p != ori_mark; ) {
            for ( auto const t: T )
            {
                if ( p == t )
                {
                    return this->make_slice( h, p );
                } 
            }
            if ( this->attribution( p ) == 1 ) { 
                auto e = this->ordering_pair( p );															
                if ( this->slice( e ) == m ) { 													
                    e = this->make_slice( h, p );																	
                    return e;
                }
            }
        }
        for ( auto p = ori_mark; --p != h; ) { 
            for ( auto const t: T )
            {
                if ( p == t )
                {
                    return this->make_slice( h, p );
                } 
            }
            if ( this->attribution( p ) == 1 ) { 
                auto e = this->ordering_pair( p );												
                if ( this->slice( e ) == m ) { 													
                    e = this->make_slice( h, p );																	
                    return e;
                }
            }
        }
    }

	// ### find_entrances_to_diff_slice			
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	find_entrances_to_diff_slice( typename base::topo_descriptor h ) 
	{
        auto Vtd = std::vector< typename base::topo_descriptor >();

        auto m = this->slice( h );

        for ( auto p = m->begin(), ori_mark = m->end(); p != ori_mark; ++p ) 
        {
            if ( this->attribution( p ) != 2 && p != h ) 
            {
                auto e = this->ordering_pair( p );															
                if ( this->slice( e ) != m ) 
                { 													
                    Vtd.push_back( p );
                }
            }
        }

        return Vtd;
    }

	// ### anticlockwise_route_through			
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	anticlockwise_route_through( typename base::topo_descriptor h, typename base::topo_descriptor t ) 
	{
        while ( true )
        {
            auto m = this->slice( h );      // assume h and t are in the same slice
            auto ori_mark = m->end();

            for ( auto p = std::next( h ); p != h; ++p ) 
            {
                if ( p == ori_mark )    { continue; }
                if ( p == t )           { return this->make_slice( h, p ); }

                if ( this->attribution( p ) == 0 ) 
                { 
                    auto e = this->ordering_pair( p );															
                    if ( this->slice( e ) == m ) 
                    { 													
                        h = this->make_slice( h, p );																	
                        break;
                    }
                }
            }
        }
    }
    
	// ### clockwise_route_through			
	
	template < typename CircularFrameSolution >
	auto circular_frame< CircularFrameSolution >::
	clockwise_route_through( typename base::topo_descriptor h, typename base::topo_descriptor t ) 
	{
        while ( true )
        {
            auto m = this->slice( h );      // assume h and t are in the same slice
            auto ori_mark = m->end();

            for ( auto p = std::prev( h ); p != h; --p ) 
            {
                if ( p == ori_mark )    { continue; }
                if ( p == t )           { return this->make_slice( h, p ); }

                if ( this->attribution( p ) == 1 ) 
                { 
                    auto e = this->ordering_pair( p );															
                    if ( this->slice( e ) == m ) 
                    { 													
                        h = this->make_slice( h, p );																	
                        break;
                    }
                }
            }
        }
    }


    // # comparasion

    template < typename CircularFrameNet, typename GeometricalPath > 
    struct comparision
    {
        CircularFrameNet    n1;
        CircularFrameNet    n2;
        GeometricalPath     gp1;
        GeometricalPath     gp2;

        comparision( CircularFrameNet const& _n1, GeometricalPath const& _gp1, CircularFrameNet const& _n2, GeometricalPath const& _gp2 ):
            n1( _n1 ), n2( _n2 ), gp1( _gp1 ), gp2( _gp2 )
        {}
       
        auto write_to_xlsx( std::string file_name ); 
        auto write_to_strm( std::ostream& os );
        auto sort_by_net( CircularFrameNet& n );
    };

    // ## write_to_xlsx 

    template < typename CircularFrameNet, typename GeometricalPath > 
    auto comparision< CircularFrameNet, GeometricalPath >::write_to_xlsx( std::string file_name )
    {
        auto workbook_name = (file_name + ".xlsx").c_str();

        lxw_workbook  *workbook  = workbook_new( workbook_name );
        lxw_worksheet *worksheet = workbook_add_worksheet( workbook, nullptr );

        // write header

        worksheet_write_string( worksheet, 0, 0, "net"          , nullptr );
        worksheet_write_string( worksheet, 0, 1, "gp1"          , nullptr );
        worksheet_write_string( worksheet, 0, 2, "gp2"          , nullptr );
        worksheet_write_string( worksheet, 0, 3, "normalization", nullptr );

        // write net and gp1

        auto Idx1 = sort_by_net( n1 ); 
        for ( int row = 0; row != n1.size(); ++row )
        {
            auto i      = Idx1[row];
            auto net    = fmt::format( "{{{} {}}}", n1[i][0], n1[i][1] ).c_str();
            auto length = boost::geometry::length( gp1[ i ] );
            worksheet_write_string( worksheet, row + 1, 0, net      , nullptr );
            worksheet_write_number( worksheet, row + 1, 1, length   , nullptr );
        }

        // write gp2

        auto Idx2 = sort_by_net( n2 ); 
        for ( int row = 0; row != n2.size(); ++row )
        {
            auto i      = Idx2[row];
            auto net    = fmt::format( "{{{} {}}}", n2[i][0], n2[i][1] ).c_str();
            auto length = boost::geometry::length( gp2[ i ] );
            worksheet_write_number( worksheet, row + 1, 2, length, nullptr );
        }

        // write normalization 
        
        for ( int row = 0; row != n1.size(); ++row )
        {
            auto f = fmt::format( "=C{0} / B{0}", row + 2 ).c_str();
            worksheet_write_formula( worksheet, row + 1, 3, f, NULL);
        }

        // write average 

        auto g1 = fmt::format( "=AVERAGE(B{}：B{})", 2, n1.size() + 1 ).c_str();
        auto g2 = fmt::format( "=AVERAGE(C{}：C{})", 2, n1.size() + 1 ).c_str();
        auto g3 = fmt::format( "=C{0} / B{0}", n1.size() + 3 ).c_str();
        worksheet_write_formula( worksheet, n1.size() + 2, 1, g1, NULL);
        worksheet_write_formula( worksheet, n1.size() + 2, 2, g2, NULL);
        worksheet_write_formula( worksheet, n1.size() + 2, 3, g3, NULL);
    
        return workbook_close(workbook);
    }

    // ## write_to_strm 

    template < typename CircularFrameNet, typename GeometricalPath > 
    auto comparision< CircularFrameNet, GeometricalPath >::write_to_strm( std::ostream& os )
    {
    }

    // ## sort_by_net 

    template < typename CircularFrameNet, typename GeometricalPath > 
    auto comparision< CircularFrameNet, GeometricalPath >::sort_by_net( CircularFrameNet& n )
    {
        for ( auto& net: n ) { if ( net[0] > net[1] ) std::swap( net[0], net[1] ); }

        auto IDX = std::vector< unsigned >( n.size() );
        std::iota( IDX.begin(), IDX.end(), 0 );
        std::sort( IDX.begin(), IDX.end(), 
            [&]( auto i, auto j )
            {
                if      ( n[ i ][ 0 ] < n[ j ][ 0 ] ) 
                {
                    return true;
                }
                else if ( n[ i ][ 0 ] == n[ j ][ 0 ] )
                {
                    return n[ i ][ 1 ] < n[ j ][ 1 ];
                }
                else 
                {
                    return false;
                }
            } 
       );
       return IDX;
    }

}
// test function declartions
void topology_test_1( void );
void topology_test_2( void );
void topology_test_3( void );
void topology_test_4( void );
void topology_test_5( void );
void topology_test_6( void );
void topology_test_6_2( void );
void topology_test_6_3( void );
void topology_test_7( void );
void topology_test_8( std::string ViaResults_FileName, std::string IoDrc_FileName );
void topology_test_9( void );
void topology_test_10( void );
void topology_test_11( std::string ViaResults_FileName, std::string IoDrc_FileName );

#endif

#ifndef _SYC_CIRCULARFRAME_H
#define _SYC_CIRCULARFRAME_H

#include <vector>
#include <tuple>
#include <list>

namespace syc::topology::v_04
{
    namespace detail
    {
        struct  topological_vertex;
    }

    // declaration
    namespace decl   
    {
        using   mapping_list            = std::list< std::vector< detail::topological_vertex > >;
        using   slice_base              = std::list< std::tuple< std::size_t, mapping_list::iterator, unsigned > >;
        using   topology_descriptor     = slice_base::iterator;
        using   slice_descriptor        = std::list< slice_base >::iterator;
        using   circular_frame_type     = std::list< std::tuple< slice_descriptor, topology_descriptor, slice_descriptor > >;
        using   circular_descriptor     = circular_frame_type::iterator;
    } 
    
    namespace detail
    {
        struct topological_vertex
        {
            std::size_t             
                id;             
            decl::circular_descriptor     
                position;       
        };

        struct slice : public decl::slice_base
        {
            public:
                decl::topology_descriptor   
                    start;
                decl::topology_descriptor   
                    end;

            private:
                decl::slice_descriptor				
                    mother;
                std::list< decl::slice_descriptor >	
                    children;

            template < typename VertexListGraph > 
                friend class slice_list;
        };
    }
   
    namespace model
    {
        template < typename VertexListGraph >
        class slice_list 
        {
            // data members
            private:
                std::vector< decl::mapping_list > 
                    _embedding_frame;

                decl::circular_frame_type
                    _circular_frame;

                std::list< detail::slice >
                    _topology_plane;

            // constructor 
            public:
                slice_list( VertexListGraph const& g );

            // topology_paths
            public:
                auto topology_paths( void );

            // add_slice 
            public:
                auto add_slice( decl::topology_descriptor const& s, decl::topology_descriptor const& t );

            // remove_slice 
            public:
                auto free_slice( decl::topology_descriptor const& h );

        };
    } 

    // # planar_directed_forset

	template < typename CoordinateType, typename PointType, typename VoronoiDiagramType >
	struct planar_directed_forset_solution :
        public std::vector< std::size_t >   // predecessor map
    {
		// ## data members
        
		private:
			std::size_t			                    _num_roots;
			std::vector<PointType>					_points;
			std::shared_ptr<VoronoiDiagramType>		_vd;

        // ## constructor

        template < typename RootContainer, typename LeafContainer, typename NetContainer >
            planar_directed_forset_solution( RootContainer&& roots, LeafContainer&& leafs, NetContainer&& nets );

        // ## num_roots and num_vertices
        
        public: 
            std::size_t num_roots( void ) { return _num_roots; }
            std::size_t num_vertices( void ) { return this->size(); }

    };

}

#endif

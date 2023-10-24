#ifndef _GNUGRAPH_H
#define _GNUGRAPH_H

#include <fstream>
#include <string>

namespace syc 
{
	namespace gnuplot
	{
		namespace v_01 
		{
			// gnugraph
			class gnugraph : public std::pair<std::ofstream, std::ofstream> // the first one is the setting script, the second one is the plot script
			{
				// copy control
				public: 
					gnugraph( std::string const& script_name ): _scriptName( script_name )
					{
						this->first.open( _scriptName + ".stg" );
						this->second.open( _scriptName + ".plt" );
						this->second << "plot 1/0 title ''\\";
					}
				// setting
				public:
					template <typename StringType>
						void set_terminal( StringType const& terminal_name );
					template <typename StringType>
						void set_output( StringType const& output_name );
					template <typename StringType>
						void set_title( StringType const& name );
					template <typename T1, typename T2>
						void set_xrange( T1 const l, T2 const u, std::string const& precision );
					template <typename T1, typename T2>
						void set_yrange( T1 const l, T2 const u, std::string const& precision );
					template <typename StringType>
						unsigned set_style_line( StringType const& styles );
					template< typename T1, typename T2, typename T3, typename T4>
						unsigned set_object_rectangle( T1 x, T2 y, T3 w, T4 h, std::string const& other_object_properties, std::string const& precision );
					template< typename T1, typename T2, typename T3>
						unsigned set_object_circle( T1 x, T2 y, T3 s, std::string const& other_object_properties, std::string const& precision );
				// plot
				public:
					template <typename StringType1, typename StringType2>
						std::shared_ptr<std::ofstream> plot_datafile_with_points( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index );
					template <typename StringType1, typename StringType2>
						std::shared_ptr<std::ofstream> plot_datafile_with_lines( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index );
					template <typename StringType1, typename StringType2>
						std::shared_ptr<std::ofstream> plot_datafile_with_labels_point( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index );

				// output
				public:
					void plot( std::string const& option ) 
					{
						this->first.close();
						this->second.close();
						std::string cmd = std::string("gnuplot ") + option + " " + _scriptName + ".stg " + _scriptName + ".plt";
						system( cmd.c_str() ); 
					}
				// clean
				public:
					void clean_all( void ) 
					{
						clean_scripts();
						clean_data_files();
					}
					void clean_scripts( void ) 
					{ 
						std::string cmd = std::string("rm ") + _scriptName + ".stg " + _scriptName + ".plt";
						system( cmd.c_str() ); 
					}
					void clean_data_files( void ) 
					{ 
						std::string cmd = "rm";
						for ( auto const& i: _data_file_names ) {
							cmd += " ";
							cmd += i;
						}
						system( cmd.c_str() ); 
					}
				// data member
				private:
					unsigned					_lineStyleSize = 0;
					unsigned					_objectSize = 0;
					std::string					_scriptName;
					std::vector<std::string>	_data_file_names;
			};
			// gnugraph setting
			template <typename StringType>
				inline void gnugraph::
				set_terminal( StringType const& terminal_name )
				{
					this->first << "set terminal " << terminal_name << "\n";
				}
			template <typename StringType>
				inline void gnugraph::
				set_output( StringType const& output_name )
				{
					this->first << "set output " << output_name<< "\n";
				}
			template <typename StringType>
				inline void gnugraph::
				set_title( StringType const& name ) 
				{
					this->first << "set title '" << name << "'\n";
				}
			template <typename T1, typename T2>
				inline void gnugraph::
				set_xrange( T1 const l, T2 const u, std::string const& precision ) 
				{
					this->first << "set xrange [" << l << precision << ":" << u << precision << "]\n";
				}
			template <typename T1, typename T2>
				inline void gnugraph::
				set_yrange( T1 const l, T2 const u, std::string const& precision ) 
				{
					this->first << "set yrange [" << l << precision <<":"<< u << precision  <<"]\n";
				}
			template <typename StringType>
				inline unsigned gnugraph::
				set_style_line( StringType const& styles )
				{
					this->first << "set style line " << ++_lineStyleSize << " " << styles << "\n";
					return _lineStyleSize;
				}
			template< typename T1, typename T2, typename T3, typename T4>
				inline unsigned gnugraph::
				set_object_rectangle( T1 x, T2 y, T3 w, T4 h, std::string const& other_object_properties, std::string const& precision ) 
				{
					this->first << "set object " << ++_objectSize << " rectangle at " << x << precision  << "," << y << precision  << " size " << w << precision  << "," << h << precision  << " "  
						  << other_object_properties << "\n";
					return _objectSize;
				}
			template< typename T1, typename T2, typename T3>
				inline unsigned gnugraph::
				set_object_circle( T1 x, T2 y, T3 s, std::string const& other_object_properties, std::string const& precision ) 
				{
					this->first << "set object " << ++_objectSize << " circle at " << x << precision  << "," << y << precision  << " size " << s << precision  << " " 
						  << other_object_properties << "\n";
					return _objectSize;
				}
			// gnugraph plot
			template <typename StringType1, typename StringType2>
				inline std::shared_ptr<std::ofstream> gnugraph::
				plot_datafile_with_points( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index )
				{
					this->second << "\n, '" << data_file_name << "' title '" << title_name << "' with points ls " << line_stlye_index << "\\";
					_data_file_names.push_back( data_file_name );
					std::shared_ptr<std::ofstream> data_file( new std::ofstream( data_file_name ) );
					return data_file;
				}	
			template <typename StringType1, typename StringType2>
				inline std::shared_ptr<std::ofstream> gnugraph::
				plot_datafile_with_labels_point( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index )
				{
					this->second << "\n, '" << data_file_name << "' title '" << title_name << "' with labels point ls " << line_stlye_index << " offset char 0.5,0.5 \\";
					_data_file_names.push_back( data_file_name );
					std::shared_ptr<std::ofstream> data_file( new std::ofstream( data_file_name ) );
					return data_file;
				}	
			template <typename StringType1, typename StringType2>
				inline std::shared_ptr<std::ofstream> gnugraph::
				plot_datafile_with_lines( StringType1 const& data_file_name, StringType2 const& title_name, unsigned line_stlye_index )
				{
					this->second << "\n, '" << data_file_name << "' title '" << title_name << "' with lines ls " << line_stlye_index << "\\";
					_data_file_names.push_back( data_file_name );
					std::shared_ptr<std::ofstream> data_file( new std::ofstream( data_file_name ) );
					return data_file;
				}
		}
		namespace v_02 
		{
		}
	}
}


#endif

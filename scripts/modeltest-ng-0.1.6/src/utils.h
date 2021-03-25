/*
  Copyright (C) 2016 Diego Darriba

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Diego Darriba <Diego.Darriba@h-its.org>,
  Heidelberg Institute for Theoretical Studies,
  Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#ifndef UTILS_H
#define UTILS_H

#include "msa.h"
#include "tree.h"
#include "global_defs.h"

#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>

#define BYTE_TO_GB 1073741824

namespace modeltest {

#ifdef WIN32
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};
#else
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '/';
    }
};
#endif

#define  ASCII_SIZE                128
#define  EOS                       0x00000200

#define  SYM_CR                    1 << 0
#define  SYM_LF                    1 << 1
#define  SYM_LFCR                  1 << 2
#define  SYM_DIGIT                 1 << 3
#define  SYM_CHAR                  1 << 4
#define  SYM_SPACE                 1 << 5
#define  SYM_TAB                   1 << 6
#define  SYM_EOF                   1 << 7
#define  SYM_UNKNOWN               1 << 8
#define  SYM_DOT                   1 << 9
#define  SYM_COLON                 1 << 10
#define  SYM_OPAREN                1 << 11
#define  SYM_CPAREN                1 << 12
#define  SYM_COMMA                 1 << 13
#define  SYM_SEMICOLON             1 << 14
#define  SYM_EQUAL                 1 << 15
#define  SYM_DASH                  1 << 16
#define  SYM_SLASH                 1 << 17
#define  SYM_PLUS                  1 << 18
#define  SYM_OBRACKET              1 << 19
#define  SYM_CBRACKET              1 << 20

#define  TOKEN_NUMBER              1 << 0
#define  TOKEN_STRING              1 << 1
#define  TOKEN_EOF                 1 << 2
#define  TOKEN_WHITESPACE          1 << 3
#define  TOKEN_NEWLINE             1 << 4
#define  TOKEN_UNKNOWN             1 << 5
#define  TOKEN_COLON               1 << 6
#define  TOKEN_OPAREN              1 << 7
#define  TOKEN_CPAREN              1 << 8
#define  TOKEN_FLOAT               1 << 9
#define  TOKEN_COMMA               1 << 10
#define  TOKEN_SEMICOLON           1 << 11
#define  TOKEN_EQUAL               1 << 12
#define  TOKEN_DASH                1 << 13
#define  TOKEN_SLASH               1 << 14
#define  TOKEN_OBRACKET            1 << 15
#define  TOKEN_CBRACKET            1 << 16

#define CONSUME(x)         while (token.tokenType & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

typedef struct
 {
   int 	        tokenType;
   const char * lexeme;
   long         len;
 } lexToken;

class Utils
{
public:
    Utils();

    /**
     * @brief Get the set of parameters supported by a tool (RAxML, PhyML, ...)
     * @param[in] tool the tool
     * @param[in] datatype the data type (dna, amino acid)
     * @return the set of supported parameters
     */
    static mt_mask_t get_parameters_from_template(template_models_t tool,
                                                  data_type_t datatype);

    /**
     * @brief Get the protein empirical matrices supported by a tool (RAxML, PhyML, ...)
     * @param[in] tool the tool
     * @param[out] n_matrices the number of matrices supported
     * @return the array of supported empirical matrices
     */
    static const mt_index_t *get_prot_matrices_from_template(template_models_t tool,
                                                             mt_size_t *n_matrices);

    /**
     * @brief Get the total number of models out of the number of matrices and set of model parameters
     * @param[in] n_matrices the number of matrices
     * @param[in] model_params the set of model parameters
     * @return the total number of models that can be constructed with the given parameters
     */
    static mt_size_t number_of_models(mt_size_t n_matrices,
                                      mt_mask_t model_params);

    /**
     * @brief Get the DNA substitution schemes supported by a tool (RAxML, PhyML, ...)
     * @param[in] tool the tool
     * @return the largest substitution scheme supported
     */
    static dna_subst_schemes_t get_dna_matrices_from_template(template_models_t tool);

    /**
     * @brief Gets the base file of a path+filename
     * @param[in] filename The filename
     * @return The name of the file, excluding the path
     */
    static std::string getBaseName(std::string const& filename);

    /**
     * @brief Estimates the required memory for P-matrices, CLVs and  scalers
     * @param[in] n_taxa The number of taxa
     * @param[in] n_sites The sequence length
     * @param[in] n_categories The number of gamma rate categories
     * @param[in] n_states The number of states (e.g., 4 nt, 20 aa)
     * @return The estimated required memory, in bytes
     */
    static size_t mem_size(unsigned int n_taxa,
                           unsigned int n_sites,
                           unsigned int n_categories,
                           unsigned int n_states);

    /**
     * @brief Allocates an uninitialized chunk of memory and checks the return value
     * @param[in] n    The number of elements
     * @param[in] size The size of each element
     */
    static void * allocate(mt_size_t n, mt_size_t size);

    /**
     * @brief Allocates a chunk of memory initialized to zero and checks the return value
     * @param[in] n    The number of elements
     * @param[in] size The size of each element
     */
    static void * c_allocate(mt_size_t n, mt_size_t size);

    /**
     * @brief Exits modeltest and prints an error message
     * @param[in] message The error message
     * @param[in] ...     Formatting arguments
     */
    static void exit_with_error(const char * message, ...) __attribute__ ((noreturn));

    /**
     * @brief Parses a file containing a set of partitions
     * @param[in] filename The file to parse
     * @return The set of partitions sorted by starting site
     */
    static partitioning_scheme_t * parse_partitions_file (std::string filename);

    /**
     * @brief Sorts a partitioning scheme according to the starting sites
     * @param[in,out] scheme The scheme to sort
     */
    static void sort_partitioning_scheme(partitioning_scheme_t & scheme);

    /**
     * @brief Count the number of set bits
     * @param[in] value the value
     * @return the number of set bits
     */
    static mt_size_t count_bits( uint32_t value);

    /**
     * @brief Parse a size variable from a string
     * @param[in] str the string representation
     * @return the size representation
     */
    static mt_size_t parse_size(const char *str);

    /**
     * @brief Parse an index variable from a string
     * @param[in] str the string representation
     * @return the index representation
     */
    static mt_index_t parse_index(const char *str);

    /**
     * @brief Converts a time (seconds) to string
     * @param[in] seconds the time in seconds
     * @return The string representation
     */
    static std::string format_time(time_t seconds);

    static std::string int_array_to_string(const int array[], int length);

    /* File utils */

    /**
     * @brief Check if a file exists
     * @param filename the file to check
     * @return true, if file with name `filename` exists
     */
    static bool file_exists(const std::string & filename);

    /**
     * @brief Check if a file is writable
     * Warning: The content of the file is removed
     * @param filename the file to check
     * @return true, if file with name `filename` is writable
     */
    static bool file_writable(const std::string & filename);

    static std::ofstream * open_file_for_writing(const std::string & filename);

    /**
     * @brief Appends a text to a file
     * Open the file, append the `text`, close the file
     * @param filename the file to write to
     * @param text the text to append
     * @return true, if text was written
     */
    static bool append_to_file(const std::string & filename,
                               const std::string & text);

    /**
     * @brief Appends a text to a file
     * @param outfile the output file stream
     * @param text the text to append
     * @return true, if text was written
     */
    static bool append_to_file(std::ofstream & outfile,
                               const std::string & text);

    /* System utils */

    /**
     * @brief Count the number of physical CPU cores
     * @return the number of physical CPU cores
     */
    static mt_size_t count_physical_cores( void );

    /**
     * @brief Count the number of logical CPU cores
     * @return the number of logical CPU cores
     */
    static mt_size_t count_logical_cores( void );

    /**
     * @brief Get the available memory
     * @return the amount of memory
     */
    static unsigned long get_memtotal( void );
};

}

#endif // UTILS_H

#pragma once
#ifndef MPI_UTILITIES_HPP
#define MPI_UTILITIES_HPP


#include <iostream>


int getSize();

int getRank();

// void MPI_print(std::ostream& stream, const std::string& message) ;
// inline void MPIsafe_print(std::ostream& stream, const std::string& message);
void MPIsafe_print(std::ostream& stream, const std::string& message);

void MPIsafe_exit(int exit_code);

void MPIsafe_barrier();
// inline void MPIsafe_barrier();

void MPIsafe_bcast_string(std::string& s, int root);

#endif
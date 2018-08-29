/******************************************************************************
 *
 * @section MODIFICATION
 *
 * Modification by Ruida Zhong for the R package VICmodel on Sep. 3th, 2018:
 * This file is largely modified to meet the requirements of CRAN and running
 * in the R runtime environment. The file is mainly modified as:
 * The macro of `log_info()`, `log_warn()`, `log_err()` are modified to adapt
 * the R environment, since `exit()` would crash the entire R process.
 * `FILE *LOG_DEST;` is added with "extern" for global variables defined in
 * head file is inavailable for C++, and definition of this variable has been
 * removed to global.cpp in /src.
 *
 * @section DESCRIPTION
 *
 * Logging macros
 *
 * @section LICENSE
 *
 * Copyright (c) 2010, Zed A. Shaw and Mongrel2 Project Contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * Neither the name of the Mongrel2 Project, Zed A. Shaw, nor the names
 *       of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#ifndef vic_log_h
#define vic_log_h

extern FILE *LOG_DEST;

#ifndef LOG_LVL
#define LOG_LVL 25
#endif

void finalize_logging(void);
void get_logname(const char *path, int id, char *filename);
void initialize_log(void);
void print_trace(void);
void setup_logging(int id, char log_path[], FILE **logfile);

// Print and error functions of R
void Rprintf(const char *, ...);
void Rf_error(const char *, ...);
void Rf_warning(const char *, ...);

// Debug Level
#if LOG_LVL < 10
#define debug Rprintf
#else
#define debug(M, ...)

#endif

// Info Level
#if LOG_LVL < 20
#define log_info Rprintf

#else
#define log_info(M, ...)
#endif

// Warn Level
#if LOG_LVL < 30
#define log_warn Rf_warning

#else
#define log_warn(M, ...)
#endif

#define log_err Rf_error

#define check_alloc_status(A, M) if (A == NULL) {Rf_error(M "%s\n", "");}

#endif

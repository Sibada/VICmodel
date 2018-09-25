#include <vic_R.h>

void Date2dmy(Date date, dmy_struct* dmy) {
  dmy->year = date.getYear();
  dmy->month = date.getMonth();
  dmy->day = date.getDay();
  dmy->dayseconds = 0;
}

unsigned short int my_str_to_agg_type(String aggstr)
  {
      if (aggstr == "avg") {
        return AGG_TYPE_AVG;
      }
      else if (aggstr == "begin") {
        return AGG_TYPE_BEG;
      }
      else if (aggstr == "end") {
        return AGG_TYPE_END;
      }
      else if (aggstr == "max") {
        return AGG_TYPE_MAX;
      }
      else if (aggstr == "min") {
        return AGG_TYPE_MIN;
      }
      else if (aggstr == "sum") {
        return AGG_TYPE_SUM;
      }
      else if (aggstr == "default") {
        return AGG_TYPE_DEFAULT;
      }
      else {
        log_err("Unknown aggregation type found: %s", aggstr.get_cstring());
      }
  }


unsigned short int my_str_to_freq(String freq)
{
  if (freq == "never") {
    return FREQ_NEVER;
  }
  else if (freq == "step") {
    return FREQ_NSTEPS;
  }
  else if (freq == "second") {
    return FREQ_NSECONDS;
  }
  else if (freq == "minute") {
    return FREQ_NMINUTES;
  }
  else if (freq == "hour") {
    return FREQ_NHOURS;
  }
  else if (freq == "day") {
    return FREQ_NDAYS;
  }
  else if (freq == "month") {
    return FREQ_NMONTHS;
  }
  else if (freq == "year") {
    return FREQ_NYEARS;
  }
  else if (freq == "date") {
    return FREQ_DATE;
  }
  else if (freq == "end") {
    return FREQ_END;
  }
  else {
    log_err("Unknown frequency flag found: %s", freq.get_cstring());
  }

}


void make_output_info(List output_infos, stream_struct **streams,
                      dmy_struct     *dmy_current)
{
  extern option_struct       options;

  size_t                     sn;
  char                       varname[MAXSTRING];
  char                       format[MAXSTRING];
  int                        type;
  double                     mult;
  unsigned short int         freq;
  int                        freq_n;
  dmy_struct                 freq_dmy;
  unsigned short int         agg_type;
  int                        n_streams;

  List                       output_info;
  CharacterVector            varnames;
  CharacterVector            aggtypes;
  int                        n_aggtypes;

  strcpy(format, OUT_ASCII_FORMAT_DEFAULT);
  mult = OUT_MULT_DEFAULT;
  type = OUT_TYPE_DEFAULT;
  n_streams = output_infos.length();
  options.Noutstreams = n_streams;

  // Allocate streams
  *streams = (stream_struct*)calloc(options.Noutstreams, sizeof(*(*streams)));
  check_alloc_status(*streams, "Memory allocation error.");

    for(sn = 0; sn < options.Noutstreams; sn++){
      output_info = as<List>(output_infos[sn]);
      varnames = as<CharacterVector>(output_info["outvars"]);

      setup_stream(&(*streams)[sn], varnames.length(), 1);
      freq = my_str_to_freq(as<String>(output_info["timescale"]));
      if (freq == FREQ_DATE) {
        if (! is<Date>(output_info["aggpar"])) {
          log_err("AGGFREQ was set to DATE but no date string was found");
        }
        Date2dmy(as<Date>(output_info["aggpar"]), &freq_dmy);
        set_alarm(dmy_current, freq, &freq_dmy, (&(*streams)[sn].agg_alarm));
      } else {
        if (! is<NumericVector>(output_info["aggpar"])) {
          freq_n = 1;
        }
        else {
          freq_n = as<int>(output_info["aggpar"]);
        }
        set_alarm(dmy_current, freq, &freq_n, (&(*streams)[sn].agg_alarm));
      }

      if(!is<CharacterVector>(output_info["aggtypes"])) {
        n_aggtypes = 0;
      } else {
        aggtypes = as<CharacterVector>(output_info["aggtypes"]);
        n_aggtypes = aggtypes.length();
      }

      for(int i = 0; i < varnames.length(); i++) {
        if(i < n_aggtypes) {
          agg_type = my_str_to_agg_type(aggtypes[i]);
        } else {
          agg_type = AGG_TYPE_DEFAULT;
        }
        strcpy(varname, ((String)varnames[i]).get_cstring());
        set_output_var(&(*streams)[sn], varname, i,
                       format, type, mult, agg_type);
      }
    }
  for (sn = 0; sn < options.Noutstreams; sn++) {
    // Allocate memory for the stream aggdata arrays
    alloc_aggdata(&(*streams)[sn]);
  }
}


int str_to_varid(String str) {
  extern metadata_struct     out_metadata[N_OUTVAR_TYPES];
  for (int varid = 0; varid < N_OUTVAR_TYPES; varid++) {
    if (str == out_metadata[varid].varname) {
      return varid;
    }
  }
  return -1;
}


List make_output_tables(List output_infos) {
  extern global_param_struct global_param;
  extern metadata_struct     out_metadata[N_OUTVAR_TYPES];

  int ncols;
  int nrows;
  int varid;

  List                       output_tables;
  List                       output_info;
  CharacterVector            varnames;
  CharacterVector            aggtypes;
  IntegerVector              varids;
  String                     tablename;
  String                     varname;
  String                     hdrname;
  int                        n_output_tables;

  n_output_tables = output_infos.length();

  for(int i = 0; i < n_output_tables; i++) {
    tablename = (String) ((CharacterVector)output_infos.names())[i];
    ncols = 0;
    output_info = as<List>(output_infos[i]);
    varnames = as<CharacterVector>(output_info["outvars"]);

    nrows = as<int>(output_info["nrow"]);
    if(nrows <= 0) nrows = global_param.nrecs;

    varids = IntegerVector(varnames.length(), 0);

    for(int j = 0; j < varnames.length(); j++) {
      varname = varnames[j];
      varid = str_to_varid(varname);
      if(varid == -1)
        log_err("set_output_var: \"%s\" was not found in the list of "
                  "supported output variable names.  Please use the exact name "
                  "listed in https://vic.readthedocs.io/en/latest/"
                  "Documentation/OutputVarList/",
                  varname.get_cstring());

      ncols += out_metadata[varid].nelem;
      varids[j] = varid;
    }

    CharacterVector header(ncols);
    int a = 0;
    for(int j = 0; j < varnames.length(); j++) {
      varid = varids[j];
      for (int k = 0; k < out_metadata[varid].nelem; k++) {
        hdrname = out_metadata[varid].varname;
        if (out_metadata[varid].nelem > 1)
          hdrname = String(sprintf<MAXSTRING>(
            "%s_%d", hdrname.get_cstring(), k+1));

        header[a] = hdrname;
        a++;
      }
    }

    NumericMatrix table(nrows, ncols);
    table.attr("header") = header;

    CharacterVector time(nrows, "");
    table.attr("time") = time;

    output_tables.push_back(table, tablename);
  }

  return output_tables;
}



void write_data(stream_struct **streams, dmy_struct     *dmy,
                List &output_tables, IntegerVector &write_row) {
  extern option_struct options;
  extern metadata_struct     out_metadata[N_OUTVAR_TYPES];

  size_t sn;
  int    col;
  int    varid;
  NumericMatrix   data_table;

  int hour, minute, second;

  for (sn = 0; sn < options.Noutstreams; sn++) {
    col = 0;
    data_table = as<NumericMatrix>(output_tables[sn]);
    if (raise_alarm(&(*streams)[sn].agg_alarm, dmy)) {
      for (int v = 0; v < (*streams)[sn].nvars; v++) {
        varid = (*streams)[sn].varid[v];

        for (int e = 0; e < out_metadata[varid].nelem; e++) {
          data_table(write_row[sn], col) =
            (*streams)[sn].aggdata[0][v][e][0];
          col ++;
        }
      }

      // Write time
      CharacterVector time = as<CharacterVector>(data_table.attr("time"));

      second = (*streams)[sn].time_bounds[0].dayseconds;
      hour = second / 3600;
      second -= hour*3600;
      minute = second / 60;
      second -= minute*60;
      time[write_row[sn]] = sprintf<MAXSTRING>("%04d-%02d-%02d %02d:%02d:%02d",
              (*streams)[sn].time_bounds[0].year,
              (*streams)[sn].time_bounds[0].month,
              (*streams)[sn].time_bounds[0].day,
              hour, minute, second);

      write_row[sn] ++;
      reset_stream((&(*streams)[sn]), dmy);
    }

  }
}

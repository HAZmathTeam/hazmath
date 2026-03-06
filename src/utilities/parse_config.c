/*
 * Configuration File Parser for AMR
 *
 * Parses files with format: variable_name{variable_value}
 * or variable_name{value value value ...} for arrays of doubles.
 * Values within braces can span multiple lines.
 * Comments start with % and continue to end of line.
 */

#include "hazmath.h"

/*
 * strdup0 - Duplicate string (safe version)
 */
static inline char* strdup0(const char* s)
{
  if (s == NULL) return NULL;
  size_t len = strlen(s) + 1;
  char* p = (char*)malloc(len);
  if (p) memcpy(p, s, len);
  return p;
}

/*
 * is_integer - Check if string represents an integer
 */
static INT is_integer(const char* str)
{
  if (str == NULL || *str == '\0') return 0;
  if (*str == '-' || *str == '+') str++;
  if (*str == '\0') return 0;
  while (*str) {
    if (!isdigit((unsigned char)*str)) return 0;
    str++;
  }
  return 1;
}

/*
 * is_double - Check if string represents a floating point number
 */
static INT is_double(const char* str)
{
  if (str == NULL || *str == '\0') return 0;
  const char* p = str;
  INT has_dot = 0, has_exp = 0, has_digits = 0;

  if (*p == '-' || *p == '+') p++;
  while (*p) {
    if (isdigit((unsigned char)*p)) {
      has_digits = 1;
    } else if (*p == '.') {
      if (has_dot || has_exp) return 0;
      has_dot = 1;
    } else if (*p == 'e' || *p == 'E') {
      if (has_exp || !has_digits) return 0;
      has_exp = 1;
      has_digits = 0;
      if (*(p + 1) == '+' || *(p + 1) == '-') p++;
    } else {
      return 0;
    }
    p++;
  }
  return has_digits;
}

/*
 * trim - Remove leading and trailing whitespace
 */
static char* trim(char* str)
{
  if (str == NULL) return NULL;
  while (isspace((unsigned char)*str)) str++;
  if (*str == '\0') return str;
  char* end = str + strlen(str) - 1;
  while (end > str && isspace((unsigned char)*end)) end--;
  end[1] = '\0';
  return str;
}

/*
 * strip_comment - Remove everything after comment char (%, #, !) to end of line
 */
static void strip_comment(char* line)
{
  char* earliest = NULL;
  char* p;

  /* Find earliest occurrence of any comment character */
  p = strchr(line, '%');
  if (p != NULL && (earliest == NULL || p < earliest)) earliest = p;

  p = strchr(line, '#');
  if (p != NULL && (earliest == NULL || p < earliest)) earliest = p;

  p = strchr(line, '!');
  if (p != NULL && (earliest == NULL || p < earliest)) earliest = p;

  if (earliest != NULL) {
    *earliest = '\0';
  }
}

/*
 * count_tokens - Count space-separated tokens in a string
 */
static INT count_tokens(const char* str)
{
  INT count = 0;
  INT in_token = 0;

  while (*str) {
    if (isspace((unsigned char)*str)) {
      in_token = 0;
    } else {
      if (!in_token) {
        count++;
        in_token = 1;
      }
    }
    str++;
  }
  return count;
}

/*
 * all_tokens_numeric - Check if all space-separated tokens are valid numbers
 */
static INT all_tokens_numeric(const char* str)
{
  char* copy = strdup0(str);
  if (copy == NULL) return 0;

  INT all_numeric = 1;
  char* token = strtok(copy, " \t\n\r,;");
  while (token != NULL && all_numeric) {
    char* endptr;
    strtod(token, &endptr);
    /* Check if entire token was consumed */
    if (*endptr != '\0') {
      all_numeric = 0;
    }
    token = strtok(NULL, " \t\n\r,;");
  }

  free(copy);
  return all_numeric;
}

/*
 * parse_double_array - Parse space-separated doubles from string
 * Returns number of doubles parsed, stores them in array
 */
static INT parse_double_array(const char* str, REAL* array, INT max_count)
{
  INT count = 0;
  char* copy = strdup0(str);
  if (copy == NULL) return 0;

  char* token = strtok(copy, " \t\n\r,;");
  while (token != NULL && count < max_count) {
    char* trimmed = trim(token);
    if (strlen(trimmed) > 0) {
      array[count++] = strtod(trimmed, NULL);
    }
    token = strtok(NULL, " \t\n\r,;");
  }

  free(copy);
  return count;
}

/*
 * get_input - Parse configuration file with name{value} format
 *
 * Format:
 *   variable_name{variable_value}   % comment ignored
 *   another_var{123}
 *   array_var{1.0 2.0 3.0
 *             4.0 5.0 6.0}
 *
 * Returns config_z structure with all parsed fields.
 * Arrays are always parsed as doubles.
 */
config_z get_input(const char* filename)
{
  config_z config;
  config.num_fields = 0;

  /* Initialize all field pointers to NULL */
  for (INT i = 0; i < MAX_FIELDS; i++) {
    config.fields[i].dbl_array = NULL;
    config.fields[i].array_len = 0;
  }

  FILE* fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "***Error***: Cannot open %s\n", filename);
    return config;
  }

  /* Read entire file into buffer, stripping comments */
  char* buffer = (char*)malloc(1024 * 1024 * 128); /* 128MB buffer */
  if (buffer == NULL) {
    fprintf(stderr, "***Error***: Cannot allocate buffer\n");
    fclose(fp);
    return config;
  }

  size_t buf_pos = 0;
  char line[MAX_STR_LEN];

  while (fgets(line, sizeof(line), fp) != NULL) {
    strip_comment(line);
    size_t len = strlen(line);
    if (buf_pos + len < 1024 * 1024 - 1) {
      memcpy(buffer + buf_pos, line, len);
      buf_pos += len;
    }
  }
  buffer[buf_pos] = '\0';
  fclose(fp);

  /* Parse name{value} pairs from buffer */
  char* ptr = buffer;

  while (*ptr && config.num_fields < MAX_FIELDS) {
    /* Skip whitespace */
    while (*ptr && isspace((unsigned char)*ptr)) ptr++;
    if (*ptr == '\0') break;

    /* Find opening brace */
    char* open_brace = strchr(ptr, '{');
    if (open_brace == NULL) break;

    /* Find matching closing brace */
    char* close_brace = strchr(open_brace, '}');
    if (close_brace == NULL) break;

    /* Extract name (handle optional = before {) */
    char name[MAX_NAME_LEN];
    size_t name_len = open_brace - ptr;

    /* Remove trailing = and whitespace from name */
    while (name_len > 0 && (ptr[name_len - 1] == '=' ||
                            isspace((unsigned char)ptr[name_len - 1]))) {
      name_len--;
    }

    if (name_len == 0 || name_len >= MAX_NAME_LEN) {
      ptr = close_brace + 1;
      continue;
    }

    strncpy(name, ptr, name_len);
    name[name_len] = '\0';
    char* trimmed_name = trim(name);

    /* Extract value (content between braces) */
    size_t value_len = close_brace - open_brace - 1;
    char* value = (char*)malloc(value_len + 1);
    if (value == NULL) {
      ptr = close_brace + 1;
      continue;
    }

    strncpy(value, open_brace + 1, value_len);
    value[value_len] = '\0';
    char* trimmed_value = trim(value);

    /* Store field */
    field_z* f = &config.fields[config.num_fields];
    snprintf(f->name, MAX_NAME_LEN, "%s", trimmed_name);
    f->dbl_array = NULL;
    f->array_len = 0;

    /* Count tokens to determine if it's an array */
    INT num_tokens = count_tokens(trimmed_value);

    if (num_tokens == 0) {
      /* Empty value - store as empty string */
      f->type = TYPE_STRING;
      f->value.str_val[0] = '\0';
    } else if (num_tokens == 1) {
      /* Single value - determine type */
      if (is_integer(trimmed_value)) {
        f->type = TYPE_INTEGER;
        f->value.int_val = strtol(trimmed_value, NULL, 10);
      } else if (is_double(trimmed_value)) {
        f->type = TYPE_DOUBLE;
        f->value.dbl_val = strtod(trimmed_value, NULL);
      } else {
        f->type = TYPE_STRING;
        strncpy(f->value.str_val, trimmed_value, MAX_STR_LEN - 1);
        f->value.str_val[MAX_STR_LEN - 1] = '\0';
      }
    } else if (all_tokens_numeric(trimmed_value)) {
      /* Multiple numeric tokens - parse as REAL array */
      f->type = TYPE_DOUBLE_ARRAY;
      f->dbl_array = (REAL*)malloc(num_tokens * sizeof(REAL));
      if (f->dbl_array != NULL) {
        f->array_len =
            parse_double_array(trimmed_value, f->dbl_array, num_tokens);
      }
    } else {
      /* Multiple non-numeric tokens - store as string */
      f->type = TYPE_STRING;
      strncpy(f->value.str_val, trimmed_value, MAX_STR_LEN - 1);
      f->value.str_val[MAX_STR_LEN - 1] = '\0';
    }

    free(value);
    config.num_fields++;
    ptr = close_brace + 1;
  }
  free(buffer);

  /* Convert all field names to lowercase */
  for (INT i = 0; i < config.num_fields; i++) {
    for (char* p = config.fields[i].name; *p; p++) {
      *p = tolower((unsigned char)*p);
    }
  }

  return config;
}
/*---------------------------------------------------------------------------*/
/*
 * config2vars_amr - Convert config_z to input_grid
 */
void config2vars_amr(config_z* config, input_grid* gz)
{
  /* Initialize with defaults (strdup so input_grid_free can free them) */
  gz->title = strdup("");
  gz->fgrid = strdup("output/3d_c.msh");
  gz->fvtu = strdup("output/3d_c.vtu");
  gz->dim = 3;
  gz->print_level = 0;
  gz->ncsys = 1;
  gz->nv = 8;
  gz->ne = 3;
  gz->nel = 1;
  gz->nf = 6;
  gz->nref = 1;
  gz->ref_type = 20;
  gz->mark_type = 0;
  gz->err_stop = -1.e-10;
  gz->num_refine_points = 1;

  /* Initialize pointers to NULL */
  gz->syslabels = NULL;
  gz->systypes = NULL;
  gz->ox = NULL;
  gz->xv = NULL;
  gz->csysv = NULL;
  gz->labelsv = NULL;
  gz->bcodesv = NULL;
  gz->xe = NULL;
  gz->seg = NULL;
  gz->mnodes = NULL;
  gz->mfaces = NULL;
  gz->data_refine_points = NULL;

  /* Temporary pointers for data arrays */
  REAL *data_coordsystems = NULL;
  REAL *data_vertices = NULL;
  REAL *data_edges = NULL;
  REAL *data_macroelements = NULL;
  REAL *data_macrofaces = NULL;
  REAL *data_refine_points = NULL;

  /* Grab values from config */
  for (INT i = 0; i < config->num_fields; i++) {
    field_z* f = &config->fields[i];
    if (!strcmp(f->name, "title")) {
      free(gz->title); gz->title = strdup(f->value.str_val);
    } else if (!strcmp(f->name, "file_grid")) {
      free(gz->fgrid); gz->fgrid = strdup(f->value.str_val);
    } else if (!strcmp(f->name, "file_vtu")) {
      free(gz->fvtu); gz->fvtu = strdup(f->value.str_val);
    } else if (!strcmp(f->name, "dimension")) {
      gz->dim = f->value.int_val;
    } else if (!strcmp(f->name, "num_coordsystems")) {
      gz->ncsys = f->value.int_val;
    } else if (!strcmp(f->name, "num_vertices")) {
      gz->nv = f->value.int_val;
    } else if (!strcmp(f->name, "num_edges")) {
      gz->ne = f->value.int_val;
    } else if (!strcmp(f->name, "num_macroelements")) {
      gz->nel = f->value.int_val;
    } else if (!strcmp(f->name, "num_macrofaces")) {
      gz->nf = f->value.int_val;
    } else if (!strcmp(f->name, "num_refinements")) {
      gz->nref = f->value.int_val;
    } else if (!strcmp(f->name, "refinement_type")) {
      gz->ref_type = f->value.int_val;
    } else if (!strcmp(f->name, "amr_marking_type")) {
      gz->mark_type = f->value.int_val;
    } else if (!strcmp(f->name, "print_level")) {
      gz->print_level = (SHORT)f->value.int_val;
    } else if (!strcmp(f->name, "num_refine_points")) {
      gz->num_refine_points = f->value.int_val;
    } else if (!strcmp(f->name, "err_stop_refinement")) {
      gz->err_stop = f->value.dbl_val;
    } else if (!strcmp(f->name, "data_coordsystems")) {
      data_coordsystems = f->dbl_array;
    } else if (!strcmp(f->name, "data_vertices")) {
      data_vertices = f->dbl_array;
    } else if (!strcmp(f->name, "data_edges")) {
      data_edges = f->dbl_array;
    } else if (!strcmp(f->name, "data_macroelements")) {
      data_macroelements = f->dbl_array;
    } else if (!strcmp(f->name, "data_macrofaces")) {
      data_macrofaces = f->dbl_array;
    } else if (!strcmp(f->name, "data_refine_points")) {
      data_refine_points = f->dbl_array;
    }
  }

  INT dim = gz->dim;
  INT i, j;

  /* data_coordsystems: [label][ox[0:dim-1]][type], row_width = dim+2 */
  if (gz->ncsys > 0 && data_coordsystems != NULL) {
    INT row_w = dim + 2;
    gz->syslabels = (INT*)malloc(gz->ncsys * sizeof(INT));
    gz->systypes = (INT*)malloc(gz->ncsys * sizeof(INT));
    gz->ox = (REAL*)malloc(gz->ncsys * dim * sizeof(REAL));
    for (i = 0; i < gz->ncsys; i++) {
      gz->syslabels[i] = (INT)data_coordsystems[i * row_w + 0];
      for (j = 0; j < dim; j++)
        gz->ox[i * dim + j] = data_coordsystems[i * row_w + 1 + j];
      gz->systypes[i] = (INT)data_coordsystems[i * row_w + dim + 1];
    }
  }

  /* data_vertices: [label][csys][xv[0:dim-1]], row_width = dim+2 */
  if (gz->nv > 0 && data_vertices != NULL) {
    INT row_w = dim + 2;
    gz->labelsv = (INT*)malloc(gz->nv * sizeof(INT));
    gz->csysv = (INT*)malloc(gz->nv * sizeof(INT));
    gz->xv = (REAL*)malloc(gz->nv * dim * sizeof(REAL));
    gz->bcodesv = (INT*)calloc(gz->nv, sizeof(INT));
    for (i = 0; i < gz->nv; i++) {
      gz->labelsv[i] = (INT)data_vertices[i * row_w + 0];
      gz->csysv[i] = (INT)data_vertices[i * row_w + 1];
      for (j = 0; j < dim; j++)
        gz->xv[i * dim + j] = data_vertices[i * row_w + 2 + j];
    }
  }

  /* data_edges: [v1][v2][divisions], row_width = 3 */
  if (gz->ne > 0 && data_edges != NULL) {
    gz->seg = (INT*)malloc(gz->ne * 3 * sizeof(INT));
    gz->xe = (REAL*)malloc(gz->ne * dim * sizeof(REAL));
    for (i = 0; i < gz->ne; i++) {
      INT v1 = (INT)data_edges[i * 3 + 0];
      INT v2 = (INT)data_edges[i * 3 + 1];
      INT div = (INT)data_edges[i * 3 + 2];
      /* ensure seg[i][0] < seg[i][1] */
      if (v1 < v2) {
        gz->seg[i * 3 + 0] = v1;
        gz->seg[i * 3 + 1] = v2;
      } else {
        gz->seg[i * 3 + 0] = v2;
        gz->seg[i * 3 + 1] = v1;
      }
      gz->seg[i * 3 + 2] = div;
    }
    /* Compute xe as midpoints of edges */
    if (gz->xv != NULL) {
      for (i = 0; i < gz->ne; i++) {
        INT v1 = gz->seg[i * 3 + 0];
        INT v2 = gz->seg[i * 3 + 1];
        for (j = 0; j < dim; j++)
          gz->xe[i * dim + j] = 0.5 * (gz->xv[v1 * dim + j] + gz->xv[v2 * dim + j]);
      }
    }
  } else {
    /* Default edge when ne=0 or data_edges is empty */
    if (gz->xv == NULL || gz->nv < 2) {
      fprintf(stderr, "***Error***: nv < 2, cannot create default edge\n");
      input_grid_free_arrays(gz);
      exit(1);
    }
    gz->ne = 1;
    gz->seg = (INT*)malloc(3 * sizeof(INT));
    gz->seg[0] = 0;
    gz->seg[1] = 1;
    gz->seg[2] = 1;
    gz->xe = (REAL*)malloc(dim * sizeof(REAL));
    /* Compute midpoint of default edge (vertices 0 and 1) */
    for (i = 0; i < dim; i++)
      gz->xe[i] = 0.5 * (gz->xv[0 * dim + i] + gz->xv[1 * dim + i]);
  }

  /* data_macroelements: row_width = 2^dim + 1 (cube vertices + material code) */
  if (gz->nel > 0 && data_macroelements != NULL) {
    INT nvcube = (1 << dim);
    INT row_w = nvcube + 1;
    gz->mnodes = (INT*)malloc(gz->nel * row_w * sizeof(INT));
    for (i = 0; i < gz->nel * row_w; i++)
      gz->mnodes[i] = (INT)data_macroelements[i];
  }

  /* data_macrofaces: row_width = 2^(dim-1) + 1 (face vertices + boundary code) */
  if (gz->nf > 0 && data_macrofaces != NULL) {
    INT nvface = (1 << (dim - 1));
    INT row_w = nvface + 1;
    gz->mfaces = (INT*)malloc(gz->nf * row_w * sizeof(INT));
    for (i = 0; i < gz->nf * row_w; i++)
      gz->mfaces[i] = (INT)data_macrofaces[i];
  }

  /* data_refine_points: row_width = dim+2 */
  if (gz->num_refine_points > 0 && data_refine_points != NULL) {
    INT row_w = dim + 2;
    INT len = gz->num_refine_points * row_w;
    gz->data_refine_points = (REAL*)malloc(len * sizeof(REAL));
    for (i = 0; i < len; i++)
      gz->data_refine_points[i] = data_refine_points[i];
  }

  /* Print input_grid elements */
  printf("\n%%input_grid params:\n"
         "  title       = '%s';\n"
         "  dim         = %6d;\n"
         "  print_level = %6d;\n"
         "  fgrid       = '%s';\n"
         "  fvtu        = '%s';\n"
         "  ncsys       = %6d;\n"
         "  nv          = %6d;\n"
         "  ne          = %6d;\n"
         "  nel         = %6d;\n"
         "  nf          = %6d;\n"
         "  nref        = %6d;\n"
         "  ref_type    = %6d;\n"
         "  mark_type   = %6d;\n"
         "  err_stop    = %9.4e;\n"
         "  num_refine_points = %d;\n",
         gz->title, gz->dim, gz->print_level, gz->fgrid, gz->fvtu,
         gz->ncsys, gz->nv, gz->ne, gz->nel, gz->nf,
         gz->nref, gz->ref_type, gz->mark_type, gz->err_stop,
         gz->num_refine_points);

  return;
}
/*-----------------------------------------------------------------------------------------------*/
/*
 * input_grid_alloc - Allocate and initialize input_grid structure
 */
input_grid* input_grid_alloc(void)
{
  input_grid* g = (input_grid*)malloc(sizeof(input_grid));
  if (g == NULL) return NULL;

  /* Initialize scalar fields with defaults */
  g->title = NULL;
  g->fgrid = NULL;
  g->fvtu = NULL;
  g->dim = 0;
  g->print_level = 0;
  g->ncsys = 0;
  g->nv = 0;
  g->ne = 0;
  g->nel = 0;
  g->nf = 0;
  g->nref = 0;
  g->ref_type = 0;
  g->mark_type = 0;
  g->err_stop = 0.0;
  g->num_refine_points = 0;

  /* Initialize all pointers to NULL */
  g->syslabels = NULL;
  g->systypes = NULL;
  g->ox = NULL;
  g->xv = NULL;
  g->csysv = NULL;
  g->labelsv = NULL;
  g->bcodesv = NULL;
  g->xe = NULL;
  g->seg = NULL;
  g->mnodes = NULL;
  g->mfaces = NULL;
  g->data_refine_points = NULL;

  return g;
}
/*-----------------------------------------------------------------------------------------------*/
/*
 * input_grid_free_arrays - Free arrays inside input_grid (for stack-allocated structs)
 */
void input_grid_free_arrays(input_grid* g)
{
  if (g == NULL) return;

  if (g->syslabels) { free(g->syslabels); g->syslabels = NULL; }
  if (g->systypes) { free(g->systypes); g->systypes = NULL; }
  if (g->ox) { free(g->ox); g->ox = NULL; }
  if (g->xv) { free(g->xv); g->xv = NULL; }
  if (g->csysv) { free(g->csysv); g->csysv = NULL; }
  if (g->labelsv) { free(g->labelsv); g->labelsv = NULL; }
  if (g->bcodesv) { free(g->bcodesv); g->bcodesv = NULL; }
  if (g->xe) { free(g->xe); g->xe = NULL; }
  if (g->seg) { free(g->seg); g->seg = NULL; }
  if (g->mnodes) { free(g->mnodes); g->mnodes = NULL; }
  if (g->mfaces) { free(g->mfaces); g->mfaces = NULL; }
  if (g->data_refine_points) { free(g->data_refine_points); g->data_refine_points = NULL; }
}
/*-----------------------------------------------------------------------------------------------*/
/*
 * free_config - Free memory allocated for arrays in config
 */
void free_config(config_z* config)
{
  if (config == NULL) return;

  for (INT i = 0; i < config->num_fields; i++) {
    if (config->fields[i].dbl_array != NULL) {
      free(config->fields[i].dbl_array);
      config->fields[i].dbl_array = NULL;
    }
    config->fields[i].array_len = 0;
  }
  config->num_fields = 0;
}

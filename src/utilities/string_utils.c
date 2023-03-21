/*! \file src/utilities/string_utils.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/17.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note: 20230321
 */
/********************************************************************/
#include "hazmath.h"

/*****************************************************************************************/
/*!
 * \fn void trim_str(char **str, INT nstrings)
 *
 * \brief function to print the ascii codes separated by | of the chars in a string s. 
 *
 * \param *s the string;
 */
void print_ascii(const char *s)
{
  size_t l=strlen(s)+1;
  const char *s_end=s+l;
  fprintf(stdout,"\nl=%ld\n|",l);
  for(; s < s_end ;++s)
    fprintf(stdout,"%03d|",*s);
  fprintf(stdout,"\n");
  return;
}
/*****************************************************************************************/
/*!
 * \fn void trim_str(char **str, INT nstrings)
 *
 * \brief function to trim all space-like chars from the beginning and the end of a string. 
 *        The new string is realloc-ed and stored at str[k], k=0:nstrings-1
 *
 * \param **str pointer to array of strings. 
 */
void trim_str(char **str, INT nstrings)
{
  const char *back,*front;
  size_t size_s,l0;
  INT k;
  for(k=0;k<nstrings;k++){
    front=str[k];
    l0=strlen(front)+1;//include the '\0' at the end
    // Trim spaces in front:
    while(isspace((unsigned char)*front)){
      front++;
    }
    if(*front == 0)  {
      // this is an empty string
      *str[k] = '\0';
      continue;
    }
    // Trim back
    back = front + strlen(front) - 1;
    while((back > front) && isspace((unsigned char)*back)) {
      back--;
    }
    back++; // ad one more at the end
    //  Set size of the trimmed string:
    size_t j,l1=l0-1,n_s=(size_t )(back-front);
    if(n_s < l1){
      size_s = n_s;
      // Copy the string at the beginning, realloc and add '\0'.
      for(j=0;j<size_s;++j){
	*(str[k]+j)=*front;
	front++;
      }
      // no need to copy anything;
    } else {
      size_s=l1;
    }
    str[k]=realloc(str[k],(size_s+1)*sizeof(char));
    *(str[k]+size_s) = '\0';
  }
  return;
}
/*****************************************************************************************/
/*!
 * \fn void dos2unix(char **str, INT nstrings)
 *
 * \brief converts \r\n to \n and single \r to \n and overwittes the input. 
 *
 * \param **str pointer to array of strings to work on. 
 */
void dos2unix(char **str, INT nstrings)
{
  const char *s0;
  char *s;
  size_t l0,size_s,j;
  INT k;
  for(k=0;k<nstrings;k++){
    l0=strlen(str[k]);
    s0=str[k];
    s=str[k];
    if(((INT )l0-1)>0){
      j=0;
      while(j<(l0-1)){
	if(*(s0+j)=='\r' && *(s0+j+1)=='\n'){
	  j+=2;
	  *s='\n';
	  s++;
	} else {
	  *s=*(s0+j);
	  j++;
	  s++;
	}
      }
      size_s=s-str[k];
      str[k]=realloc(str[k],(size_s+1)*sizeof(char));
      *(str[k]+size_s)='\0';
      s=str[k];
      while(*s!='\0'){
	if(*s=='\r') *s='\n';
	s++;
      }
    }
  }
  return;
}

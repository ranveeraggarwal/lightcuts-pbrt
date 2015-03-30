/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_CORE_PBRTPARSE_HPP_INCLUDED
# define YY_YY_CORE_PBRTPARSE_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    STRING = 258,
    ID = 259,
    NUM = 260,
    LBRACK = 261,
    RBRACK = 262,
    ACCELERATOR = 263,
    ACTIVETRANSFORM = 264,
    ALL = 265,
    AREALIGHTSOURCE = 266,
    ATTRIBUTEBEGIN = 267,
    ATTRIBUTEEND = 268,
    CAMERA = 269,
    CONCATTRANSFORM = 270,
    COORDINATESYSTEM = 271,
    COORDSYSTRANSFORM = 272,
    ENDTIME = 273,
    FILM = 274,
    IDENTITY = 275,
    INCLUDE = 276,
    LIGHTSOURCE = 277,
    LOOKAT = 278,
    MAKENAMEDMATERIAL = 279,
    MATERIAL = 280,
    NAMEDMATERIAL = 281,
    OBJECTBEGIN = 282,
    OBJECTEND = 283,
    OBJECTINSTANCE = 284,
    PIXELFILTER = 285,
    RENDERER = 286,
    REVERSEORIENTATION = 287,
    ROTATE = 288,
    SAMPLER = 289,
    SCALE = 290,
    SHAPE = 291,
    STARTTIME = 292,
    SURFACEINTEGRATOR = 293,
    TEXTURE = 294,
    TRANSFORMBEGIN = 295,
    TRANSFORMEND = 296,
    TRANSFORMTIMES = 297,
    TRANSFORM = 298,
    TRANSLATE = 299,
    VOLUME = 300,
    VOLUMEINTEGRATOR = 301,
    WORLDBEGIN = 302,
    WORLDEND = 303,
    HIGH_PRECEDENCE = 304
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 157 "core/pbrtparse.yy" /* yacc.c:1909  */

char string[1024];
float num;
ParamArray *ribarray;

#line 110 "core/pbrtparse.hpp" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_CORE_PBRTPARSE_HPP_INCLUDED  */

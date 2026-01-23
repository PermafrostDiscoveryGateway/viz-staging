{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 HelveticaNeue-Bold;\f1\fnil\fcharset0 .SFNS-Regular;\f2\fnil\fcharset0 .AppleSystemUIFontMonospaced-Regular;
\f3\fnil\fcharset0 .SFNS-Bold;\f4\froman\fcharset0 TimesNewRomanPSMT;}
{\colortbl;\red255\green255\blue255;\red14\green14\blue14;}
{\*\expandedcolortbl;;\cssrgb\c6700\c6700\c6700;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\sl324\slmult1\pardirnatural\partightenfactor0

\f0\b\fs28 \cf2 Inputs
\f1\b0 \cf2 \
\pard\tqr\tx100\tx260\li260\fi-260\sl324\slmult1\sb240\partightenfactor0
\cf2 	\'95	
\f2 \cf2 input_path
\f1 \cf2 : Path to input vector dataset (.shp or .gpkg)\
	\'95	
\f2 \cf2 output_path
\f1 \cf2 : Path to output GeoPackage\
	\'95	h3_
\f2 \cf2 res
\f1 \cf2 : H3 resolution (1\'9615)\
	\'95	
\f2 \cf2 land_polygons_path
\f1 \cf2  (optional): Land/coastline polygon dataset path\
	\'95	
\f2 \cf2 area_epsg
\f1 \cf2 : Equal-area EPSG code for area computations (default 
\f2 \cf2 6933
\f1 \cf2 )\
\pard\tqr\tx100\tx260\li260\fi-260\sl324\slmult1\sb240\partightenfactor0
\cf2 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\sl324\slmult1\pardirnatural\partightenfactor0

\f0\b \cf2 Output
\f1\b0 \
\pard\tqr\tx100\tx260\li260\fi-260\sl324\slmult1\sb240\partightenfactor0
\cf2 	\'95	
\f2 \cf2 out_gdf
\f1 \cf2 : H3 grid (polygons) with per-cell summary attributes written to 
\f2 \cf2 output_path
\f1 \cf2  (GeoPackage)\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\sl324\slmult1\pardirnatural\partightenfactor0

\f3\b\fs30 \cf2 Algorithm 1: build_h3_summary(input_path, output_path, h3_res, land_polygons_path=None, area_epsg=6933)
\f1\b0\fs28 \
\pard\tqr\tx260\tx420\li420\fi-420\sl324\slmult1\sb240\partightenfactor0

\f4 \cf2 	1.	
\f0\b Read input
\f1\b0 \
\pard\tqr\tx660\tx820\li820\fi-820\sl324\slmult1\sb240\partightenfactor0

\f4 \cf2 	1.	
\f2 gdf \uc0\u8592  gpd.read_file(input_path)
\f1 \

\f4 	2.	If 
\f2 gdf.crs
\f1  is missing: raise error\

\f4 	3.	If 
\f2 gdf.crs.to_epsg() != 4326
\f1 : 
\f2 gdf \uc0\u8592  gdf.to_crs(epsg=4326)}
/*  Editing routine for EBS crab survey data from the by-leg, by-vessel crab */
/*  specimen tables.
/*  Output file name should indicate vessel, leg, cruise year, and times run */
/*  (e.g. akk2016_leg1_crab_edit1.lst).                                      */
/*  Script written by:  Claire Armistead                                     */
/*  Date created: 5/10/2003                                                  */
/*  Last revised: 8/23/2013 to add check for small king crab w/clutch_size > 0  */
/*                and for large king crab coded 0-0-0 (immature)             */
/*  Revised by:  Claire Armistead                                            */

accept xcrab prompt 'Enter crab table name:  '
accept xyear prompt 'Enter cruise year:  '
accept xoutput prompt 'Enter output file name:  '

clear breaks
set verify off
set heading on
set pages 80
set pause off
set termout off
set feedback off

break on species_code skip 1

spool &xoutput
prompt  Output from edit crab for &xcrab, &xyear

column vessel format 999
column cruise format 999999
column species_code format 9999999
column sex format 999
column shell_condition format 99999
column egg_color format 999999999
column egg_condition format 999999999
column clutch_size format 999999
column length format 999999
column width format 999999
column chela_height format 99999
column disease_code format 99999999999999
column disease_dorsal format 99999999
column disease_ventral format 99999999
column disease_legs format 99999999
column sampling_factor format 99999.999999

/* Inventory routine  */

prompt
prompt Inventory of vessel and cruise
   select vessel, cruise, count(*) from &xcrab
   where cruise like '&xyear%'
   group by vessel, cruise;
prompt
prompt Inventory of species codes and sex
   select vessel, cruise, species_code, sex, count(*) from &xcrab
   where cruise like '&xyear%'
   group by vessel, cruise, species_code, sex
   order by vessel,cruise,species_code,sex;
prompt
prompt Inventory of shell condition
   select distinct species_code, sex, shell_condition, count(*) from &xcrab
   where cruise like '&xyear%'
   group by species_code,sex,shell_condition
   order by species_code,sex,shell_condition;
prompt
prompt Check for invalid egg color codes
   select distinct vessel,haul,species_code,sex,egg_color,
   count(*) from &xcrab
   where cruise like '&xyear%' and sex in (2,4)
   and egg_color not in (0,2,3,4,5,6,9)
   group by vessel,haul,species_code,sex,egg_color
   order by vessel,haul,species_code,sex,egg_color;
prompt
prompt  Check for invalid egg condition codes
   select distinct vessel,haul,species_code,sex,egg_condition,
   count(*) from &xcrab
   where cruise like '&xyear%' and sex in (2,4)
   and egg_condition not in (0,1,2,3,4,5,9)
   group by vessel,haul,species_code,sex,egg_condition
   order by vessel,haul,species_code,sex,egg_condition;
prompt
prompt Check for invalid clutch_size codes
   select distinct vessel,haul,species_code,sex,clutch_size,
   count(*) from &xcrab
   where cruise like '&xyear%' and sex in (2,4)
   and clutch_size = 8 or clutch_size > 9
   group by vessel,haul,species_code,sex,clutch_size
   order by vessel,haul,species_code,sex,clutch_size;
prompt

/* Check for questionable egg codes */

prompt Check for questionable egg codes for female king crab
   select distinct vessel,haul,species_code,sex,shell_condition shell,
   egg_color,egg_condition,clutch_size
   from &xcrab
   where sex = 2 and cruise like '&xyear%' and species_code in (69322,69323)
   and shell_condition = 0 and egg_condition = 1
   or sex = 2 and cruise like '&xyear%' and species_code in (69322,69323)
   and shell_condition = 1 and egg_condition > 1
   or sex = 2 and cruise like '&xyear%' and species_code = 69322
   and shell_condition in (3,4,5) and egg_condition = 1
   or sex = 2 and cruise like '&xyear%' and species_code in (69322,69323)
   and shell_condition = 1 and egg_condition >= 2
   order by vessel,haul,species_code;
prompt
prompt Check for questionable egg codes for Tanner and snow crab
   select distinct vessel,haul,species_code,sex,shell_condition shell,
   egg_color,egg_condition,clutch_size
   from &xcrab
   where sex = 2 and cruise like '&xyear%' and species_code in (68560,68580,68590)
   and shell_condition = 1 and egg_condition >= 2
   or sex = 2 and cruise like '&xyear%' and species_code in (68560,68580,68590)
   and shell_condition = 0 and egg_condition = 1 
   order by vessel,haul,species_code;

/* Inventory of shell condition and egg codes */
prompt
prompt Inventory of shell condition and egg codes for female crabs
   select distinct species_code,shell_condition,egg_color,
   egg_condition,clutch_size,count(*) from &xcrab
   where sex = 2 and cruise like '&xyear%'
   group by species_code,shell_condition,egg_color,egg_condition,clutch_size
   order by species_code,shell_condition,egg_color,egg_condition,clutch_size;

/* Check egg codes - females and males*/
prompt
prompt Any females without egg codes?
   select distinct vessel, cruise, haul, species_code,
   egg_color,egg_condition,clutch_size from &xcrab
   where (cruise like '&xyear%' and sex = 2) and (egg_color is null
   or egg_condition is null or clutch_size is null)
   order by vessel,haul,species_code;
prompt
prompt Any males with egg codes?
   select distinct vessel, cruise, haul, species_code,
   egg_color,egg_condition,clutch_size from &xcrab
   where (cruise like '&xyear%' and sex = 1) and (egg_color is not null or
   egg_condition is not null or clutch_size is not null)
   order by vessel,haul,species_code;
prompt

/* Null lengths or width */
prompt
prompt Check for null lengths for king and hair crabs
   select vessel,haul,species_code,sex,length,width
   from &xcrab
   where cruise like '&xyear%' and species_code in (69310,69322,69323,69400)
   and length is null
   order by vessel,haul,species_code;
prompt
prompt Check for null widths for Tanner crabs
   select vessel,haul,species_code,sex,length,width
   from &xcrab
   where cruise like '&xyear%' and species_code in (68560,68580,68590)
   and width is null
   order by vessel,haul,species_code;
prompt

/* Check for small female king crab with clutch_size */

prompt Check for tiny crab with clutch_size for king crabs
   select vessel,haul,species_code,sex,length,clutch_size
   from &xcrab
   where cruise like '&xyear%' and species_code in (69322,69323)
   and length < 65
   and clutch_size > 0
   order by vessel,haul,species_code;
prompt

/* Check for small female Tanner and snow crab with clutch_size */

prompt Check for tiny crab with clutch_size for Chionoecetes crabs
   select vessel,haul,species_code,sex,length,clutch_size
   from &xcrab
   where cruise like '&xyear%' and species_code in (68560,68580,68590)
   and length < 30
   and clutch_size > 0
   order by vessel,haul,species_code;
prompt
--/* Check for large female king crab coded immature */

--prompt Check for large crab coded 0-0-0 for red king crab
--   select vessel,haul,species_code,sex,length,clutch_size
--   from &xcrab
--   where cruise like '&xyear%' and species_code in (2,3)
--   and length > 100
--   and egg_color = 0
--   and egg_condition = 0
--   and clutch_size = 1
--   order by vessel,haul,species_code;
--prompt

/* Check for tiny crab with old shell condition */

prompt Check for tiny crab with old shell for Chionoecetes crabs
   select vessel,haul,species_code,sex,width,shell_condition
   from &xcrab
   where cruise like '&xyear%' and species_code in (68560,68580,68590)
   and width < 30
   and shell_condition > 2
   order by vessel,haul,species_code;
prompt
prompt Check for tiny crab with old shell for king crabs
   select vessel,haul,species_code,sex,length,shell_condition
   from &xcrab
   where cruise like '&xyear%' and species_code in (69322,69323)
   and length < 60
   and shell_condition > 2
   order by vessel,haul,species_code;
prompt   
      
/* Min and max lengths and widths */
prompt
prompt Check the minimum and maximum crab sizes by species_code.
prompt Make sure crab that should have lengths do not have widths and vice versa
   --select species_code,sex,haul,min(length),max(length),min(width),max(width) 
   --from &xcrab
   --where cruise like '&xyear%' and length <> 999 or width <> 999
   --group by species_code,sex,haul
   --order by species_code,sex,haul
 
   select species_code, sex, min(length)min_length,
   max(length)max_length, min(width)min_width, max(width)max_width
   from &xcrab
   where cruise like '&xyear%' and length <> 999 or width <> 999
   group by species_code, sex
   order by species_code,sex; 
prompt

/* chela_height height */
prompt Check the minimum and maximum chela_height height by species_code
   select species_code, length, width, min(chela_height),max(chela_height)
   from &xcrab
   where cruise like '&xyear%' and chela_height is not null
   group by species_code,length,width
   order by species_code,length,width; 
prompt

/* Any chela_height entries for females? */
prompt Check for female crab with chela_height
prompt
   select vessel,haul,species_code,chela_height
   from &xcrab
   where cruise like '&xyear%' and sex = 2
   and chela_height is not null
   order by vessel,haul,species_code;
prompt

/* Any merus_length entries for females? */
prompt Check for female crab with merus_length data
prompt 
   select vessel,haul,species_code,merus_length
   from &xcrab
   where cruise like '&xyear%' and sex = 2
   and merus_length is not null
   order by vessel,haul,species_code;
prompt

/* Min and max merus_length length */
prompt Check the minimum and maximum merus_length entries
   select species_code,sex,min(merus_length),max(merus_length)
   from &xcrab
   where cruise like '&xyear%' and merus_length is not null
   group by species_code,sex
   order by species_code,sex;
prompt

/* Pathology code check */
prompt
prompt Check the special sample codes for black mat with no % coverage
   select vessel,haul,species_code,disease_code,
   disease_dorsal,disease_ventral,disease_legs 
   from &xcrab
   where cruise like '&xyear%'
   and (disease_code = 1 and (disease_dorsal is null and disease_ventral is null
   and disease_legs is null))
   order by vessel,haul,species_code;
prompt
prompt Check the special sample codes for bitter crab with entries in % cov
   select vessel,haul,species_code,disease_code,
   disease_dorsal,disease_ventral,disease_legs
   from &xcrab
   where cruise like '&xyear%' 
   and (disease_code = 2 and (disease_dorsal is not null
   or disease_ventral is not null or disease_legs is not null))
   order by vessel,haul,species_code;
prompt
prompt Check for null special sample with entries in disease_dorsal,2,or 3
   select vessel,haul,species_code,disease_code,
   disease_dorsal,disease_ventral,disease_legs
   from &xcrab
   where cruise like '&xyear%'
   and (disease_code is null and (disease_dorsal is not null 
   or disease_ventral is not null or disease_legs is not null))
   order by vessel,haul,species_code;
prompt
prompt Check for weird special sample codes (>9)
   select vessel,haul,species_code,disease_code
   from &xcrab
   where cruise like '&xyear%'
   and disease_code > 9
   order by vessel,haul,species_code;

/* Max sample factor */
prompt
prompt Check the maximum sampling factor by species_code
   select species_code, max(sampling_factor)maxsamplefactor from &xcrab
   where cruise like '&xyear%'
   group by species_code
   order by species_code;
prompt Make sure the minimum sampling factor is 1.0
   select vessel,cruise,haul,species_code,sampling_factor from &xcrab
   where sampling_factor < 1.0000;
   
-- Check station id's
prompt
prompt Check station id
	select haul,station from &xcrab
	where station not in
	(select station from centers)
	order by station;

prompt

prompt End of edit crab routine.
spool off
set termout on

drop table check_crab_number_detail;

create table check_crab_number_detail as
select vessel,cruise,haul,species_code,round(sum(sampling_factor))number_crab_spec_table
from winter_survey2023_crab
group by vessel,cruise,haul,species_code
order by vessel,cruise,haul,species_code;

drop table check_crab_number;

create table check_crab_number as
select a.vessel,a.cruise,a.haul,a.species_code,number_caught,number_crab_spec_table
from check_crab_number_detail a, winter_survey2023_crabsummary b
where a.vessel = b.vessel
and a.cruise = b.cruise
and a.haul = b.haul
and a.species_code = b.species_code
order by a.vessel,a.cruise,a.haul,species_code;

select vessel,cruise,haul,species_code,number_caught,number_crab_spec_table from check_crab_number
where number_caught <> number_crab_spec_table;



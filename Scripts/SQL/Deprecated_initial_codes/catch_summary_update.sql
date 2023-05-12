alter table winter_survey2023_crabsummary
add NUMBER_SAMPLED NUMBER;


alter table winter_survey2023_crabsummary
add NUMBER_CAUGHT NUMBER;


update winter_survey2023_crabsummary
set number_sampled = number_crab;

update winter_survey2023_crabsummary
set number_caught = weight;                  -- This will need to be corrected for final protocol

alter table winter_survey2023_crabsummary
drop column weight;

alter table winter_survey2023_crabsummary
drop column number_crab;


commit;
-- Haul = pot #


drop table winter2023_haul;

create table winter2023_haul
 (
  HAUL                                   NUMBER,
  START_TIME                              TIMESTAMP,
  END_TIME                               TIMESTAMP,
  SOAK_TIME                              N3MBER,
  LATITUDE                         NUMBER,
  LONGITUDE                        NUMBER                          
);

COMMIT;
alter table winter2023_haul
add survey_year    NUMBER;

update winter2023_haul set
survey_year = 2023;

alter table winter2023_haul
add cruise    NUMBER;

update winter2023_haul set
cruise = 202301;


alter table winter2023_haul
add vessel   NUMBER;

update winter2023_haul set
vessel = 162;
COMMIT;

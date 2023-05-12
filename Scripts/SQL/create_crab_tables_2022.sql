drop table winter_survey2023_crab;

create table winter_survey2023_crab 
 (SPECIES_CODE       NUMBER(6,0), 
	SEX                NUMBER(2,0), 
	LENGTH             NUMBER(4,1), 
	WIDTH              NUMBER(4,1), 
	SHELL_CONDITION    NUMBER(2,0), 
	EGG_COLOR          NUMBER(2,0), 
	EGG_CONDITION      NUMBER(2,0), 
	CLUTCH_SIZE        NUMBER(2,0), 
	MERUS_LENGTH       NUMBER(4,1), 
	CHELA_HEIGHT       NUMBER(4,1), 
	DISEASE_CODE       NUMBER(2,0), 
	DISEASE_DORSAL     NUMBER(2,0), 
	DISEASE_VENTRAL    NUMBER(2,0), 
	DISEASE_LEGS       NUMBER(2,0), 
	WEIGHT             NUMBER(8,2), 
	SAMPLING_FACTOR    NUMBER(12,5) NOT NULL ENABLE, 
	VESSEL             NUMBER(4,0) NOT NULL ENABLE, 
	CRUISE             NUMBER(6,0) NOT NULL ENABLE, 
	HAUL               NUMBER(4,0) NOT NULL ENABLE, 
	STATION            VARCHAR2(10 BYTE), 
	ID                 VARCHAR2(50 BYTE), 
	COMMENTS           VARCHAR2(2000 BYTE),
  NOTES           VARCHAR2(2000 BYTE)
   ) ;
   
drop table winter_survey2023_crabsummary;

create table winter_survey2023_crabsummary
  (VESSEL NUMBER, 
	 CRUISE NUMBER, 
	 HAUL NUMBER, 
	 COMMON_NAME VARCHAR2(30 BYTE), 
	 SPECIES_CODE NUMBER, 
	 WEIGHT NUMBER, 
	 NUMBER_CRAB NUMBER
   );


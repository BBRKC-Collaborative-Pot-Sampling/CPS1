library(RCurl)
library(RSelenium)
library(rvest)

url       <-"http://forum.axishistory.com/memberlist.php"
pgsession <-html_session(url)

pgform    <-html_form(pgsession)[[2]]

filled_form <- html_form_set(pgform,
                          "username" = "username", 
                          "password" = "password")

session_submit(pgsession,filled_form)
memberlist <- jump_to(pgsession, "http://forum.axishistory.com/memberlist.php")

page <- html(memberlist)

usernames <- html_nodes(x = page, css = "#memberlist .username") 

data_usernames <- html_text(usernames, trim = TRUE) 

######
url       <-"https://sftp.afsc.noaa.gov/ThinClient/WTM/public/#/login"   ## page to spider
pgsession <-html_session(url)               ## create session
pgform    <-html_form(pgsession)[[1]]     ## pull form from session

filled_form <- set_values(pgform,
                          `ctl00$Header2$HeaderTop1$tbUsername` = "myemail@gmail.com", 
                          `ctl00$Header2$HeaderTop1$tbPassword` = "mypassword")


read_html(url) -> html

html %>%
html_element(xpath = "/body/div[1]/div/div[2]/div/div/form/div[1]")



  html_element("body") %>%
  html_nodes("div") %>%
  html_attr("class")
  
  
  html_element("ui-view")
  html_nodes("login-form ng-pristine ng-invalid-required")

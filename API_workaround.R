library(popler)
library(dplyr)

# service functions ----------------------------------------------------------

# connect to the database
popler_connector <- function(dbname=NULL, host=NULL, port=NULL, user=NULL, password=NULL, silent=TRUE) {
  
  if (!requireNamespace("RPostgreSQL", quietly = TRUE)) {
    stop("RPostgreSQL package required to connect to postgres db", call. = FALSE)
  }
  user <- if(is.null(user)){
    if(identical(Sys.getenv("TRAVIS"), "true")){"postgres"} else {""}
  } else user
  con <- RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(),
                                host     = if(is.null(host))     "" else host,
                                dbname   = if(is.null(dbname))   "" else dbname,
                                user     = user,
                                password = if(is.null(password)) "" else password,
                                port     = if(is.null(port))     "" else port)
  #info <- RPostgreSQL::dbGetInfo(con)
  #dbplyr::src_sql("postgres", con, info=info, disco=popler:::popler_disconnector(con,"postgres",silent))
  #src_sql("postgres", con, info=info, disco=popler:::popler_disconnector(con,"postgres",silent))
}

# database open function
db_open <- function(){
  popler_connector(dbname="popler_3" ,
                   host="ec2-54-214-212-101.us-west-2.compute.amazonaws.com" ,
                   port=5432 ,
                   user="other_user" ,
                   password="bigdata" ,
                   silent=TRUE)
}

# get query 
query_get <- function(connection, query){
  # accepts a connection and a string query input
  # outputs a dataframe
  return(dplyr::tbl(connection, dbplyr::sql(query)) %>% data.frame())
}


# get datasets by proj ID --------------------------------------------------------

# open connection
conn      <- db_open()

# produce efficient query (uses WITH clause)
efficienty_query <- function( proj_id ){
  
  paste0('
         WITH project_taxa_table AS (
         select authors, authors_contact,   -- 126
         datatype,
         spatial_replication_level_1_label, spatial_replication_level_2_label, spatial_replication_level_3_label, spatial_replication_level_4_label, spatial_replication_level_5_label, proj_metadata_key, structured_type_1, structured_type_2, structured_type_3, structured_type_4,
         sppcode, genus, species, taxa_table.taxa_table_key
         FROM   taxa_table            
         JOIN site_in_project_table ON taxa_table.site_in_project_taxa_key =site_in_project_table.site_in_project_key 
         JOIN project_table         ON site_in_project_table.project_table_fkey = project_table.proj_metadata_key 
         JOIN study_site_table      ON site_in_project_table.study_site_table_fkey =study_site_table.study_site_key 
         JOIN lter_table            ON study_site_table.lter_table_fkey = lter_table.lterid 
         WHERE proj_metadata_key = ',proj_id,' )
         
         SELECT authors, authors_contact, year, day, month, sppcode, genus, species, datatype, spatial_replication_level_1_label, spatial_replication_level_1, spatial_replication_level_2_label, spatial_replication_level_2, spatial_replication_level_3_label, spatial_replication_level_3, spatial_replication_level_4_label, spatial_replication_level_4, spatial_replication_level_5_label, spatial_replication_level_5, proj_metadata_key, structure_type_1, structure_type_2, structure_type_3, structure_type_4, count_table.treatment_type_1, count_table.treatment_type_2, count_table.treatment_type_3, covariates , count_observation 
         FROM project_taxa_table
         JOIN count_table ON count_table.taxa_count_fkey = project_taxa_table.taxa_table_key  
         UNION ALL 
         
         SELECT authors, authors_contact, year, day, month, sppcode, genus, species, datatype, spatial_replication_level_1_label, spatial_replication_level_1, spatial_replication_level_2_label, spatial_replication_level_2, spatial_replication_level_3_label, spatial_replication_level_3, spatial_replication_level_4_label, spatial_replication_level_4, spatial_replication_level_5_label, spatial_replication_level_5, proj_metadata_key, structure_type_1, structure_type_2, structure_type_3, structure_type_4, biomass_table.treatment_type_1, biomass_table.treatment_type_2, biomass_table.treatment_type_3, covariates , biomass_observation 
         FROM project_taxa_table
         JOIN biomass_table ON biomass_table.taxa_biomass_fkey = project_taxa_table.taxa_table_key
         UNION ALL  
         
         SELECT authors, authors_contact, year, day, month, sppcode, genus, species, datatype, spatial_replication_level_1_label, spatial_replication_level_1, spatial_replication_level_2_label, spatial_replication_level_2, spatial_replication_level_3_label, spatial_replication_level_3, spatial_replication_level_4_label, spatial_replication_level_4, spatial_replication_level_5_label, spatial_replication_level_5, proj_metadata_key, structure_type_1, structure_type_2, structure_type_3, structure_type_4, percent_cover_table.treatment_type_1, percent_cover_table.treatment_type_2, percent_cover_table.treatment_type_3, covariates , percent_cover_observation 
         FROM project_taxa_table
         JOIN percent_cover_table ON percent_cover_table.taxa_percent_cover_fkey = project_taxa_table.taxa_table_key
         UNION ALL  
         
         SELECT authors, authors_contact, year, day, month, sppcode, genus, species, datatype, spatial_replication_level_1_label, spatial_replication_level_1, spatial_replication_level_2_label, spatial_replication_level_2, spatial_replication_level_3_label, spatial_replication_level_3, spatial_replication_level_4_label, spatial_replication_level_4, spatial_replication_level_5_label, spatial_replication_level_5, proj_metadata_key, structure_type_1, structure_type_2, structure_type_3, structure_type_4, individual_table.treatment_type_1, individual_table.treatment_type_2, individual_table.treatment_type_3, covariates , individual_observation 
         FROM project_taxa_table
         JOIN individual_table ON individual_table.taxa_individual_fkey = project_taxa_table.taxa_table_key
         UNION ALL
         
         SELECT authors, authors_contact, year, day, month, sppcode, genus, species, datatype, spatial_replication_level_1_label, spatial_replication_level_1, spatial_replication_level_2_label, spatial_replication_level_2, spatial_replication_level_3_label, spatial_replication_level_3, spatial_replication_level_4_label, spatial_replication_level_4, spatial_replication_level_5_label, spatial_replication_level_5, proj_metadata_key, structure_type_1, structure_type_2, structure_type_3, structure_type_4, density_table.treatment_type_1, density_table.treatment_type_2, density_table.treatment_type_3, covariates , density_observation 
         FROM project_taxa_table
         JOIN density_table ON density_table.taxa_density_fkey = project_taxa_table.taxa_table_key')
  
}

# download the "bad" density dataset (number 12)
#init_t    <- Sys.time()
#query_get(conn, efficienty_query( 88 ))
#Sys.time() - init_t


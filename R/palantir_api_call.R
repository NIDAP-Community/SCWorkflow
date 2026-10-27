#' palantir_api_call
#' Utility function from 3D tSNE Coordinate Template from v 75#' 
#' 
#' @param service The NIDAP API service to call
#' @param path The path of NIDAP API service
#' @param token NIDAP user toekn.
#' @param data Data to be uploaded with NIDAP API calls.
#' @param method Method to be used, including POST, GET, and DELETE 
#' 
#' @importFrom httr add_headers content DELETE GET POST
#' 
#' @export
#' 
#' @return return the content of API calls
#' 
palantir_api_call <- function(service, path, token, data, method) { 
  base_url <- 'https://gypsum.palantircloud.com/'
  url <- paste(base_url, service, path, sep="")
  auth_header <- paste("Bearer", token, sep=" ")
  headers = c("Content-Type"="application/json","Authorization"=auth_header)
  request_body_json <- data
  print(url)
  print(request_body_json)
  
  # Post method
  if(method == 'POST') {
    response <- POST(url,
                     body=request_body_json,
                     encode="json",
                     add_headers(.headers=headers))
  }
  
  # Get method
  if(method == 'GET') {
    response <- GET(url,
                    encode="json",
                    add_headers(.headers=headers))
    print(response)
  }
  
  # Delete method, USE IT CAUCIOUSLY
  if(method == 'DELETE') {
    response <- DELETE(url,
                       encode="json",
                       add_headers(.headers=headers))
    print(response)
  }
  
  return(content(response,'text', encoding = "UTF-8"))
}

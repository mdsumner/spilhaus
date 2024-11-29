
#' Project points to Spilhaus
#'
#' Input a matrix of longitude latitude coordinates and return Spilhaus coordinates. 
#' 
#' Note that this transformation is not available via standard libraries, and was
#' implemented in R and Python by Ricardo T. Lemos. 
#' 
#' @param x matrix of longlat 
#'
#' @return matrix of projected coordinates
#' @export
#'
#' @examples
#' spilhaus(cbind(0, 0))
spilhaus <- function(x){

  longitude <- x[,1L, drop = TRUE]
  latitude <- x[,2L, drop = TRUE]
  
  # constants (https://github.com/OSGeo/PROJ/issues/1851)
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  # coordinates in radians
  lon = longitude * pi / 180
  lat = latitude * pi / 180

  # conformal latitude, in radians
  lat_c = -pi / 2 + 2 * atan(
    tan(pi/4 + lat/2) * ((1 - e * sin(lat)) / (1 + e * sin(lat))) ^ (e / 2)
  )

  # transformed lat and lon, in degrees
  lat_s = 180 / pi * asin(sin(alpha) * sin(lat_c) - cos(alpha) * cos(lat_c) * cos(lon - lambda_0))
  lon_s = 180 / pi * (
    beta + atan2(
      cos(lat_c) * sin(lon - lambda_0),
      (sin(alpha) * cos(lat_c) * cos(lon - lambda_0) + cos(alpha) * sin(lat_c))
    )
  )

  # projects transformed coordinates onto plane (Adams World in a Square II)
  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projected = reproj::reproj_xy(cbind(lon_s, lat_s), source = "EPSG:4326", target = adams_ws2)
  adams_x = projected[,1]
  adams_y = projected[,2]
  spilhaus_x = -(adams_x + adams_y) / sqrt(2)
  spilhaus_y = (adams_x - adams_y) / sqrt(2)

  return(cbind(x = spilhaus_x, y = spilhaus_y)) #, adams_x, adams_y, lon_s, lat_s))
}


#' Project points from Spilhaus to longlat 
#'
#' Input a matrix of Spilhaus coordinates and return longlat coordinates. 
#' 
#' Note that this transformation is not available via standard libraries, and was
#' implemented in R and Python by Ricardo T. Lemos. 
#' 
#' @param x  matrix of coordinates in Spilhaus 
#'
#' @return matrix of longlat coordinates
#' @export
#'
#' @examples
#' spilhaus_lonlat(cbind(0, 0))
spilhaus_lonlat <- function(x) {
  spilhaus_x <- x[,1L, drop = TRUE]
  spilhaus_y <- x[,2L, drop = TRUE]
  
  # constants
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  adams_x = (spilhaus_y - spilhaus_x) * sqrt(2) / 2
  adams_y = - (spilhaus_y + spilhaus_x) * sqrt(2) / 2

  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

  projected = reproj::reproj_xy(cbind(adams_x, adams_y), source = adams_ws2, target = "EPSG:4326")
  lon_s = projected[,1]
  lat_s = projected[,2]

  #transformed coords in radians
  lon_s_rad = lon_s * pi / 180
  lat_s_rad = lat_s * pi / 180

  # conformal latitude
  lat_c = asin(sin(alpha) * sin(lat_s_rad) + cos(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta))

  # longitude, in radians
  lon = lambda_0 + atan2(
    cos(lat_s_rad) * sin(lon_s_rad - beta),
    sin(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta) - cos(alpha) * sin(lat_s_rad)
  )

  # latitude (iterative formula from https://mathworld.wolfram.com/ConformalLatitude.html)
  lat = lat_c
  for (i in 0:9) {
    lat = -0.5 * pi + 2 * atan(
      tan(pi / 4 + lat_c / 2) *
        ((1 + e * sin(lat)) / (1 - e * sin(lat))) ^ (e / 2)
    )
  }

  # coordinates in degrees
  longitude = ((lon * 180 / pi + 180) %% 360) - 180
  latitude = lat * 180 / pi

  return(cbind(x = longitude, y = latitude))
}

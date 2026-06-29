# output.R
# Ergebnisse in xlsx schreiben

schreibe_output <- function(output, art = "Deckungsstock", pfad_output) {
  
  cat("Schreibe Output...\n")
  
  # Liste in Vektor umwandeln
  output_vec <- unlist(output)
  
  # Als data.frame für openxlsx
  df <- data.frame(Ergebnis = output_vec, stringsAsFactors = FALSE)
  
  # Dateiname
  if (art == "Deckungsstock") {
    dateiname <- file.path(pfad_output, "Input_Deckungsstock.xlsx")
  } else {
    dateiname <- file.path(pfad_output, "Input_Fonds.xlsx")
  }
  
  # Schreiben
  openxlsx::write.xlsx(df, file = dateiname, colNames = FALSE)
  
  cat("Output gespeichert:", dateiname, "\n")
}
setClass(
  # Class name
  "Autocart",

  # Slot type
  representation(
    name = "character",
    birth = "Date"
  ),

  # Initialize slots
  prototype = list(
    name = as.character(NULL),
    birth = as.Date(as.character(NULL))
  )
)

#' Interactive version of the uncertainty illustration
#' @export interactive.uncertainty.illustration
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return NULL
#' @examples
#' \donttest{
#' interactive.uncertainty.illustration()
#' }
interactive.uncertainty.illustration <- function(log.level = NULL) {
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  options(shiny.reactlog = TRUE)
    ## ====== Global variables ======
    rv <- shiny::reactiveValues(logger = c())

    # ---- the following needs to be reset when data changed ----
    current.fit <- NULL # store the fit so that no need to refit everytime
    current.compare.fit <- NULL # store the fit so that no need to refit everytime

    compare.method.names <- c() # store the compare.method.name so produce correct clicking effect
    fit.counter <- 0 # this global variable comes with the is.fit.button function
    use.fit.counter <- 0

    # ----log ----
    log.to.app <- function(line) {
      cat(line)
      rv$logger <- c(shiny::isolate(rv$logger), utils::capture.output(cat(line, sep = "")))
    }
    futile.logger::flog.appender(log.to.app)

    ## ===== helper functions =====

    clean.up.after.updated.dataset <- function(session) {
      shinyjs::reset("compare.method.names")
      current.fit <<- NULL
      current.compare.fit <<- NULL
    }

    is.fit.button <- function(fit.button.count) { # use to distinguish if calling of "new.fit.update.react" is from fit.button
      if (is.null(fit.button.count) || fit.button.count == 0 || fit.button.count == fit.counter) {
        return(FALSE)
      } else{
        fit.counter <<- fit.counter + 1
        return(TRUE)
      }
    }
    is.use.fit.button <- function(use.fit.button.count) { # use to distinguish if calling of "new.fit.update.react" is from fit.button
      if (is.null(use.fit.button.count) || use.fit.button.count == 0 || use.fit.button.count == use.fit.counter) {
        return(FALSE)
      } else{
        use.fit.counter <<- use.fit.counter + 1
        return(TRUE)
      }
    }

    ## ====== UI ======
    ui <- shiny::fluidPage(
      shinyjs::useShinyjs(),
      shiny::titlePanel("Uncertainty illustration"),

      # ---- data tab ----
      shiny::tabsetPanel(
        # data tab
        shiny::tabPanel("Data",
                 shiny::sidebarLayout(
                   shiny::sidebarPanel(

                     # input for the new method
                     shiny::h3("Input data"),

                     shiny::radioButtons(inputId = "dataType", "data source", choices = c("input data", "in memory data", "example data"), selected = "example data"),

                     shiny::tags$div(title = "check the box if your data has intercept",
                                     shiny::checkboxInput(inputId = "intercept", "Intercept", FALSE)
                     ),

                     shiny::conditionalPanel(condition = "input.dataType == 'input data'",
                                             shiny::tags$div(title = "choose a file that from your computer",
                                                             shiny::fileInput(inputId = "data.file", "choose File",
                                                         multiple = TRUE,
                                                         accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                      ),
                                      shiny::h3("Data format"),

                                      shiny::tags$div(title = "check the box if your data has header",
                                                      shiny::checkboxInput(inputId = "header", "Header", TRUE)
                                      ),

                                      shiny::radioButtons(inputId = "sep", "separator",
                                                   choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                                                   selected = ","),

                                      shiny::radioButtons(inputId = "quote", "quote",
                                                   choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"),
                                                   selected = '"'),

                                      shiny::textInput("na.str", "strings to be interpreted as missing values"),

                                      shiny::selectInput("response.name", "select response name",
                                                  c())
                     ),
                     shiny::conditionalPanel(condition = "input.dataType == 'in memory data'",
                                             shiny::tags$div(title = "choose an object from memory in r",
                                                             shiny::selectInput(inputId = "data.object.name", label = "object name",
                                                                                choices = c("", setdiff(ls(".GlobalEnv"), utils::lsf.str(".GlobalEnv"))), selected = "")),
                                             shiny::selectInput("object.response.name", "select response name",
                                                                                                                                                    c())
                     ),

                     shiny::conditionalPanel(condition = "input.dataType == 'example data'",
                                             shiny::h3("Example data"),
                                             shiny::selectInput(inputId = "example.name", label = "data set",
                                                  choices = c("", "mtcars", "boston"), selected = "")
                     ),
                     shiny::actionButton(inputId = "data.button",
                                  label = "use data")
                   ),


                   shiny::mainPanel(
                     DT::DTOutput(outputId = "df.contents"),
                     shiny::uiOutput(outputId = "data.text")
                   )
                 )
        ),

        # ---- fit tab ----
        shiny::tabPanel("Model fit",
                        shiny::sidebarLayout(
                          shiny::sidebarPanel(
                            shiny::radioButtons(inputId = "fitType", "fit source", choices = c("new fit", "in memory fitted result"), selected = "new fit"),

                            shiny::conditionalPanel(condition = "input.fitType == 'new fit'",
                     # input for the new method
                            shiny::h3("Fitting the new method"),

                            shiny::selectInput(inputId = "method.names",
                                 label = "variable selection methods: ",
                                 choices = c("lasso", "elastic", "relaxo", "mcp", "scad"),
                                 selected = "", multiple = T),

                            shiny::sliderInput(inputId = "B",
                                 label = "B: ",
                                 min = 5,
                                 max = 400,
                                 value = 50),

                            shiny::actionButton(inputId = "fit.button",
                                  label = "fit")),
                     shiny::conditionalPanel(condition = "input.fitType == 'in memory fitted result'",
                                             shiny::tags$div(title = "choose the fitted result from memory in r",
                                                             shiny::selectInput(inputId = "fit.object.name", label = "fitted object name",
                                                                                choices = c("", setdiff(ls(".GlobalEnv"), utils::lsf.str(".GlobalEnv"))), selected = "")),
                                             shiny::actionButton(inputId = "use.fit.button",
                                                                 label = "use")),

                     shiny::br(),
                     shiny::br(),

                     shiny::h3("Comparing with individual variable selection methods"),

                     shiny::checkboxGroupInput(inputId = "compare.method.names",
                                        label = "variable selection methods to compare: ",
                                        choices = c("lasso", "elastic", "relaxo", "mcp", "scad"),
                                        selected = "", inline = T),

                     shiny::h3("Plotting parameters"),

                     # shiny::radioButtons(inputId = "plot.type",
                     #                     label = "boxplot type:",
                     #                     choices = c("conditonal only", "with unconditional", "with violin"),
                     #                     selected = "conditonal only"),

                     shiny::checkboxGroupInput(inputId = "plot.type",
                                               label = "boxplot type:",
                                               choices = c("with unconditional", "with violin"),
                                               selected = "", inline = T),

                     shiny::sliderInput(inputId = "q",
                                 label = "q: ",
                                 min = 0,
                                 max = 20,
                                 value = 0),
                     shiny::sliderInput(inputId = "level",
                                        label = "confidence level: ",
                                        min = 50,
                                        max = 99.99,
                                        value = 90)
                   ),
                   shiny::mainPanel(
                     shiny::uiOutput(outputId = "data.source.text"),
                     shiny::plotOutput(outputId = "new.ci.plot", height = "800px"),
                     shiny::uiOutput(outputId = "fit.text"),
                     shiny::uiOutput(outputId = "compare.text")
                     # plotOutput(outputId = "oracle_ci_plot", height = "400px")
                   )
                 )
        ),
        # ---- log tab ----
        shiny::tabPanel("Log",
                        shiny::selectInput(inputId = "log.level", label = "log level", choices = c("debug", "info", "warn", "error"), selected = "info"),
                        shiny::uiOutput(outputId = "log")
        ),
        shiny::tabPanel("Help",
                        shiny::h3("Interactive tool to visualize variable selection uncertainty"),
                        shiny::br(),
                        shiny::p("This Shiny App visualizes variable selection uncertainty for a given dataset. Please follow the instruction below on how to upload and fit the data."),
                        shiny::br(),
                        shiny::tabsetPanel(
                          shiny::tabPanel("Upload data",
                                          shiny::br(),
                                          shiny::tags$ul(
                                            shiny::tags$li(
                                              shiny::p("Use your own data:"),
                                              shiny::tags$ul(
                                                shiny::tags$li('Select "input data" as your data source'),
                                                shiny::tags$li('Choose your own file by pressing the "Browse" button. If your dataset is in a right format, the first few lines of the observations will be shown on the right hand side of the window'),
                                                shiny::tags$li("Update the data format of the data:"),
                                                shiny::tags$ul(
                                                  shiny::tags$li('If your data has an header, check the "header"'),
                                                  shiny::tags$li('Select the right separator used in your data under "separator"'),
                                                  shiny::tags$li('Select the right quotation mark used in your data under "quote"'),
                                                  shiny::tags$li('Input the strings used to represent missing data in the box under "strings to be interpreted as missing values"'),
                                                  shiny::tags$li('Select the response name (Y) in the list under "select response name"')
                                  ),
                                  shiny::tags$li('Once the data format is set, click the "use data" button to upload the data for model fitting')
                                )
                              ),
                              shiny::br(),
                              shiny::tags$li(
                                shiny::p("Use the in memery data:"),
                                shiny::tags$ul(
                                  shiny::tags$li('Select "in memory data" as your data source'),
                                  shiny::tags$li("Select the in memory R object you want to use. The first few lines of the observations will be shown on the right hand side of the window"),
                                  shiny::tags$li('Select the response name (Y) in the list under "select response name"'),
                                  shiny::tags$li('Click the "use data" button to upload the data for model fitting')
                                )
                              ),
                              shiny::br(),
                              shiny::tags$li(
                                shiny::p("Use the example datasets:"),
                                shiny::tags$ul(
                                  shiny::tags$li('Select "example data" as your data source'),
                                  shiny::tags$li("Select the example data your want to use. The first few lines of the observations will be shown on the right hand side of the window"),
                                  shiny::tags$li('Click the "use data" button to upload the data for model fitting')
                                )
                              )
                            )
                   ),
                   shiny::tabPanel("Fit the data",
                            shiny::br(),
                            shiny::tags$ul(
                              shiny::tags$li("Fit the new method"),
                              shiny::tags$ul(
                                shiny::tags$li("Make sure you have already uploaded the data before attempting to fit the data"),
                                shiny::tags$li(
                                  shiny::p("Use the in memery fit:"),
                                  shiny::tags$ul(
                                    shiny::tags$li('Select "in memory fitted result" as your fit source'),
                                    shiny::tags$li("Select the in memory R csuv object you want to use"),
                                    shiny::tags$li('Once you select R object, you can click the "use" button. If the fit data is in correct type, the uncertainty illustration graph is shown on the right hand side of the window'))),
                                shiny::tags$li(
                                  shiny::p("Have new fit:"),
                                  shiny::tags$ul(
                                    shiny::tags$li('Select the variable selection methods that you want to use in the new method under "variable selection method". You can choose more than one'),
                                    shiny::tags$li('Select the number of subsampling datasets used in the new method by sliding the bar under "B"'),
                                    shiny::tags$li('Once you select the methods and B, you can click the "fit" button to fit the new method. After the fitting is finished, the uncertainty illustration graph is shown on the right hand side of the window')))

                              ),
                              shiny::br(),
                              shiny::tags$li("Compare with other variable selection methods"),
                              shiny::tags$ul(
                                shiny::tags$li("Select the variable selection methods that you want to compare with. You can choose more than one"),
                                shiny::tags$li("Once the fitting is finished, the selections of each individual variable are show on the graph as points")
                              ),
                              shiny::br(),
                              shiny::tags$li("Variations on the uncertainty illustration"),
                              shiny::tags$ul(
                                shiny::tags$li("Select the type of plot: confidence interval-like plot or boxplot"),
                                shiny::tags$li('Select "q", which represents the q-th percentile fitted model with the lowest mse used in the new method'),
                                shiny::tags$li('If confidence interval like plot is chosen, the confidence level can be varied by using the sliding bar under "confidence level"')
                              )
                            )
                   )
                 )
        ),
        shiny::tabPanel("About",
                        shiny::br(),
                        shiny::h4("About"),
                        shiny::p("This app is an interactive version of the graphical tool proposed by 'Exploiting disagreement between high-dimensional variable selectors for uncertainty visualization'. The graphical tool visualizes the variable selection uncertainty via the method introduced by the paper."),
                        shiny::p("The aim of this app is to provide users a taste of how the tool can help visualizing uncertainties, and how the tool may be helpful for their data analysis. In this app, users can try out the visualization tool using their own data and compare the new method with some popular variable selection methods. Users can also use this app to visualize the variable selection uncertainty of the data they have and have some idea which variables are possibly more important."),
                        shiny::p("Note that this shiny app is only a 'teste': it is designed to process datasets that are relatively small. This app may not be efficient or have enough compacity to handle some larger datasets. For users who want to use the method and visualization tool to analyse a big dataset, please use the csuv in the 'CSUV' R package"),
                        shiny::br(),
                        shiny::h4("Author"),
                        shiny::p("Christine Yuen (yuenl@lse.ac.uk)")
        )
      )
    )

    ## ====== server ======
    server <- function(input, output, session) {
      # ---- update buttons / tabs according to input ----
      # observeEvent(input$dataType, {
      #   reset("example.name")
      #   reset("data.file")
      #   clean.up.after.updated.dataset(session = session)
      #   })

      shiny::observe({
        df.col.names <- colnames(file.raw.data.reactive())
        shiny::updateSelectInput(session, "response.name",
                          label = paste("Select input label", length(df.col.names)),
                          choices = df.col.names,
                          selected = utils::head(df.col.names, 1)
        )
      })

      shiny::observe({
          df.col.names <- colnames(in.memory.raw.data.reactive())
          shiny::updateSelectInput(session, "object.response.name",
                                   label = paste("Select input label", length(df.col.names)),
                                   choices = df.col.names,
                                   selected = utils::head(df.col.names, 1)
          )
      })

      # ---- get raw data (for display) ----
      file.raw.data.reactive <- shiny::reactive({
        if (!is.null(input$data.file)) {
          d <- utils::read.csv(input$data.file$datapath,
                               header = input$header,
                               sep = input$sep,
                               quote = input$quote,
                               na.strings = input$na.str)
          d <- d[stats::complete.cases(d), ]
          return(d)
        }
        return(NULL)
      })

      in.memory.raw.data.reactive <- shiny::reactive({
        if (!is.null(input$data.object.name) && input$data.object.name != "") {
          d <- get(input$data.object.name)
          if (is.null(dim(d))) {
            shiny::showNotification("please select an valid object", type = "error")
            futile.logger::flog.error("no valid object selected")
          } else{
            return(d)
          }
        }
        return(NULL)
      })

      example.raw.data.reactive <- shiny::reactive({
        d <- NULL
        if (input$example.name == "mtcars") {
          d <- datasets::mtcars
        } else if (input$example.name == "boston") {
          # library(MASS)
          d <- MASS::Boston[, c(14, 1:13)]
        }
        return(d)
      })

      raw.data.reactive <- shiny::reactive({
        if (input$dataType == "example data") {
          return(example.raw.data.reactive())
        } else if (input$dataType == "in memory data") {
          return(in.memory.raw.data.reactive())
        }else{
          return(file.raw.data.reactive())
        }
      })

      # ---- get data (for fit) ----
      example.data.reactive <- shiny::reactive({
        if (is.null(input$example.name) || input$example.name == "") {
          shiny::showNotification("please select an example data set", type = "error")
          futile.logger::flog.error("no example data set selected")
          return(NULL)
        }

        return(list(x = as.matrix(example.raw.data.reactive()[, -1]),
                     y = example.raw.data.reactive()[, 1],
                     data.name = shiny::isolate(input$example.name)))
      })

      in.memory.data.reactive <- shiny::reactive({
        # check if valid dataset is selected
        if (is.null(input$data.object.name) || input$data.object.name == "") {
          shiny::showNotification("please select an object", type = "error")
          futile.logger::flog.error("no object selected")
          return(NULL)
        }
        if (is.null(dim(get(input$data.object.name)))) {
          shiny::showNotification("please select an valid object", type = "error")
          futile.logger::flog.error("no valid object selected")
          return(NULL)
        }

        # set x and y
        i <- which(names(raw.data.reactive()) == shiny::isolate(input$object.response.name))
        return(list(y = raw.data.reactive()[, i],
                     x = as.matrix(raw.data.reactive()[, -i]),
                     data.name = shiny::isolate(input$data.object.name)))
      })

      file.data.reactive <- shiny::reactive({
        # check if valid dataset is selected
        if (is.null(input$data.file)) {
          shiny::showNotification("please select a valid data set", type = "error")
          futile.logger::flog.error("no data set selected")
          return(NULL)
        }

        # set x and y
        i <- which(names(raw.data.reactive()) == shiny::isolate(input$response.name))
        return(list(y = raw.data.reactive()[, i],
                     x = as.matrix(raw.data.reactive()[, -i]),
                     data.name = input$data.file$name))
      })

      data.reactive <- shiny::eventReactive(input$data.button, {
        d <- NULL
        if (!is.null(input$data.button) && input$data.button != 0) {
          if (input$dataType == "example data") {
            d <- example.data.reactive()
          } else if (input$dataType == "in memory data") {
            d <- in.memory.data.reactive()
          } else{
            d <- file.data.reactive()
          }
          # reset
          if (!is.null(d)) {
            clean.up.after.updated.dataset(session)
            futile.logger::flog.info("uploaded the data")
            return(c(d, list(intercept = input$intercept)))
          }
        }
        return(NULL)
      }, ignoreNULL = FALSE)

      # ---- fit data (new method) ----
      new.fit.update.reactive <- shiny::eventReactive({
        input$use.fit.button
        input$fit.button
        data.reactive()}, {
          is.from.fit.button <- is.fit.button(input$fit.button)
          is.from.use.fit.button <- is.use.fit.button(input$use.fit.button)
          if (is.from.fit.button || is.from.use.fit.button) {
            # check validity of input
            if (is.null(data.reactive())) {
              shiny::showNotification("please upload data", type = "error")
              futile.logger::flog.error("no data uploaded")
              return(NULL)
            }
            # req(isolate(data.reactive()))
            if (is.from.fit.button) {
              if (shiny::isolate(is.null(input$method.names))) {
                shiny::showNotification("please select at least one method", type = "error")
                futile.logger::flog.error("no variable method selected")
                shiny::updateCheckboxGroupInput(session, inputId = "compare.method.names", selected = "")
                shiny::updateCheckboxGroupInput(session, inputId = "plot.type", selected = "")
                return(NULL)
              }
              # req(input$method.names)

              futile.logger::flog.info("start fitting the new method")
              shiny::withProgress(message = "Getting cross-validation fit", value = 0, {
                method.names <- shiny::isolate(input$method.names)
                B <- shiny::isolate(input$B)
                current.fit <<- get.csuv.unique.fit(X = shiny::isolate(data.reactive())$x,
                                                   Y = shiny::isolate(data.reactive())$y,
                                                   intercept = shiny::isolate(data.reactive())$intercept,
                                                   method.names = method.names,
                                                   fit.percent = 0.5,
                                                   B = B,
                                                   current.fit = current.fit, num.core = 1)
              })
              rv$fit.counter <- shiny::isolate(rv$fit.counter) + 1 # not used?
              futile.logger::flog.info("finish fitting the new method")
            } else{
              if (shiny::isolate(is.null(input$fit.object.name)) || shiny::isolate(input$fit.object.name == "")) {
                shiny::showNotification("please select an object", type = "error")
                futile.logger::flog.error("no all.fits object")
                return(NULL)
              }
              obj <- get(shiny::isolate(input$fit.object.name))
              if (is.csuv.fit(obj) && !is.null(obj)) {
                current.fit <<- obj$all.fits
              } else if (is.csuv.all.fits(obj)) {
                current.fit <<- obj
              } else{
                shiny::showNotification("please select an valid all.fits object (either a csuv fit or the all.fits from a csuv fit)", type = "error")
                futile.logger::flog.error("no valid all.fits object")
                return(NULL)
              }
              rv$use.fit.counter <- shiny::isolate(rv$use.fit.counter) + 1 # not used?
            }
            return(current.fit)
          } else{
            return(NULL)
          }
        }, ignoreNULL = F)

      new.fit.reactive <- shiny::reactive({
        # rv$fit.counter
        if (!is.null(new.fit.update.reactive())) {
          method.names <- NULL
          if (shiny::isolate(input$fitType == "new fit")) {
            method.names <- shiny::isolate(input$method.names)
          } else{
            method.names <- names(new.fit.update.reactive()[[1]])
          }
          return(get.csuv.final.mod(X = shiny::isolate(data.reactive()$x),
                                    Y = shiny::isolate(data.reactive()$y),
                                    intercept = shiny::isolate(data.reactive())$intercept,
                                    unique.fit = new.fit.update.reactive(),
                                    selection.criterion = "mse",
                                    coef.est.method = lm.ols,
                                    q = input$q,
                                    method.names = method.names,
                                    B = shiny::isolate(input$B)))
        }
        return(NULL)
      })

      compare.fit.update.reactive <- shiny::eventReactive(input$compare.method.names, {
        if (is.null(input$compare.method.names)) {
          compare.method.names <<- c()
          return(NULL)
        } else{
          # isolate(source("compare_method.r", local = T))
          if (is.null(data.reactive())) {
            shiny::showNotification("please upload data", type = "error")
            shiny::updateCheckboxGroupInput(session, inputId = "compare.method.names", selected = "")
            shiny::updateCheckboxGroupInput(session, inputId = "plot.type", selected = "")# seems will call the first if again
            futile.logger::flog.error("no data uploaded")
            return(NULL)
          }
          if (is.null(new.fit.reactive())) {
            shiny::showNotification("please fit the new method first", type = "error")
            shiny::updateCheckboxGroupInput(session, inputId = "compare.method.names", selected = "")
            shiny::updateCheckboxGroupInput(session, inputId = "plot.type", selected = "")# seems will call the first if again
            futile.logger::flog.error("no fitted before comparing")
            return(NULL)
          }
          i <- which(compare.method.names %in% input$compare.method.names)
          j <- which(!(input$compare.method.names %in% compare.method.names))
          compare.method.names <<- c(compare.method.names[i], input$compare.method.names[j])

          current.compare.fit <<- get.compare.fit(x = shiny::isolate(data.reactive()$x),
                                                  y = shiny::isolate(data.reactive()$y),
                                                  intercept = shiny::isolate(data.reactive()$intercept),
                                                  method.names = input$compare.method.names,
                                                  current.compare.fit = current.compare.fit)
          return(current.compare.fit)
        }

      }, ignoreNULL = FALSE)

      compare.fit.reactive <- shiny::reactive({
        if (!is.null(compare.fit.update.reactive())) {
          # validate(need(length(input$compare.method.names)>1 || input$compare.method.names!= "", "No compare method selected"))
          # compare.methods = lapply(shiny::isolate(input$compare.method.names), function(method.name) compare.fit.update.reactive()[[method.name]])
          compare.methods <- lapply(shiny::isolate(input$compare.method.names), function(method.name) compare.fit.update.reactive()[method.name, , drop = FALSE])
          names(compare.methods) <- shiny::isolate(input$compare.method.names)
          compare.methods <- do.call(rbind, compare.methods)
          return(compare.methods)
        }
        return(NULL)
      })



      # ---- data output ----
      output$df.contents <- DT::renderDT({
        return(raw.data.reactive())
      }, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))

      output$data.text <- shiny::renderUI({
        if (!is.null(data.reactive())) {
          return(shiny::HTML(paste(shiny::br(),
                                    "The current uploaded data set for fitting is",
                                    data.reactive()$data.name,
                                    "(uploaded on", Sys.time(), ")")))
        }else{
          return(NULL)
        }

      })

      # ---- plot output ----
      output$data.source.text <- shiny::renderUI({
        if (!is.null(data.reactive())) {
          return(shiny::HTML(paste(shiny::br(), "Using data set ",
                                    data.reactive()$data.name, "uploaded on",
                                    Sys.time())))
        }else{
          return(shiny::HTML("please upload a dataset"))
        }
      })

      output$new.ci.plot <- shiny::renderPlot({
        plot.choice <- input$plot.type
        suppressWarnings(csuv.plot.helper(new.fit = new.fit.reactive(),
                                          with.unconditional = ("with unconditional" %in% plot.choice),
                                          with.violin = ("with violin" %in% plot.choice),
                                          level = 1 - input$level / 100, # note the lhs is the significant level x whereas rhs is confidence level x%
                                          print.compare.method.points = TRUE,
                                          compare.method.fit = compare.fit.reactive(),
                                          compare.method.names = compare.method.names)) #
      })

      output$fit.text <- shiny::renderUI({
        if (!is.null(new.fit.reactive())) {
          return(shiny::HTML(paste("Graph last update on fitting the new method:", Sys.time())))
        } else{
          return(NULL)
        }
      })
      output$compare.text <- shiny::renderUI({
        if (!is.null(compare.fit.reactive())) {
          return(shiny::HTML(paste("Last update on comparing with other methods:", Sys.time())))
        }
      })
      # ---- log output ----
      output$log <- shiny::renderUI({
        log.levels <- c("debug", "info", "warn", "error")
        log.level <- which(log.levels == input$log.level)

        lines <- ""
        for (line in rv$logger) {
          line.log.level <- which(log.levels == tolower(sub(" .*", "", line)))
          if (line.log.level >= log.level) {
            lines <- paste(lines, shiny::br(), line)
          }
        }
        return(shiny::HTML(lines))
      })
    }
    return(shiny::shinyApp(onStart = start.func, ui = ui, server = server))
}

start.func <- function() {
  futile.logger::flog.threshold(futile.logger::INFO)
  futile.logger::flog.info("Starting the application")
  csuv.env$is.shiny <- TRUE

  shiny::onStop(function() {
    csuv.env$is.shiny <- FALSE
    futile.logger::flog.info("Stopping the application")
  })
}

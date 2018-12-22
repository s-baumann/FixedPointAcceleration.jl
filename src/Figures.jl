"""
A function for plotting the change in each vector per iterate.
### Takes
 * Inputs - The Inputs matrix returned by FixedPoint.
 * Outputs - The Outputs matrix returned by FixedPoint.
 * ConvergenceVector - The Convergence vector returned by fixedpoint.
 * ConvergenceSigFig - The number of significant figures convergence values should be shown with in the plot header.
 * ShowInputs - A boolean describing whether to show the inputs in the figure.
 * ShowOutputs - A boolean describing whether to show the outputs in the figure.
 * ShowPrevious - A boolean describing whether to show the previous inputs/outputs in the figure.
 * xaxis - A vector for meaningful values corresponding to the input/output values. By default the indexes of each input/output vector are used.
 * secondhold - If this is -1 or less then all plotting happens immediately. This means that very quickly only the last iterate will be seen so it makes sense to do it only if FromIterate and ToIterate indicate only one iterate. If this is 0 then a user can click through each figure. If this is positive then it describes how many seconds to pause between frames. By default this is 0 so a user must click through each frame.
 * FromIterate - This describes what iterate to show first.
 * ToIterate - This describes what iterate to show last.
### Returns
 * This function returns nothing. It just shows a plot in the console.
### Examples
Inputs = seq(1,10)
Function = function(x){ cos(x) }
A = FixedPoint(Function, Inputs, Algorithm = "Anderson")
ChangePerIterate(A$Inputs, A$Outputs, A$Convergence)
 # Any now to have it play one frame every half a second starting from the nineth iterate
ChangePerIterate(A$Inputs, A$Outputs, A$Convergence, secondhold = 0.5, FromIterate = 9)
"""
function ChangePerIterate(Inputs, Outputs, ConvergenceVector = c(), ConvergenceSigFig = 5, ShowInputs = TRUE, ShowOutputs = TRUE, ShowPrevious  = TRUE, xaxis = c(), secondhold = 0, FromIterate = 1, ToIterate = c()){

  if (ShowInputs + ShowOutputs < 0.5){stop("It is not possible to use this function without showing inputs or outputs (as nothing will be drawn). Set ShowInputs and/or ShowOutputs to TRUE.")}
  if (ShowInputs){ylab = "Inputs (in red)."}
  if (ShowOutputs){ylab = "Outputs (in blue)."}
  if (all(ShowInputs, ShowOutputs)){ylab = "Inputs (in red), Outputs (in blue)."}

  if (is.null(xaxis)){xaxis = 1:size(Inputs)[1]}
  if (is.null(ToIterate)){ToIterate = size(Inputs)[2]}

  ytop = max(Inputs[ , FromIterate:ToIterate],Outputs[, FromIterate:ToIterate])
  ybot = min(Inputs[ , FromIterate:ToIterate],Outputs[, FromIterate:ToIterate])

  for (i in FromIterate:ToIterate){
    if (is.null(ConvergenceVector)) { Title = paste0("Fixed Point Convergence. Iterate:", i, ".")} else {Title = paste0("Fixed Point Convergence. Iterate:", i, ". Convergence: ", NicePrint(ConvergenceVector[i],ConvergenceSigFig))}

    if (ShowPrevious){graphics::plot(c(min(xaxis), max(xaxis)), c(ytop, ybot), type = "p", col = 0, xlab = "", ylab = ylab, main = Title, sub = "The previous iterate is represented by the open circles.")
    } else {graphics::plot(c(min(xaxis), max(xaxis)), c(ytop, ybot), type = "p", col = 0, xlab = "", ylab = ylab, main = Title)}
    if (ShowInputs){graphics::points(xaxis, Inputs[,i], type = "p", pch = 19, col = "red")}
    if (ShowOutputs){graphics::points(xaxis, Outputs[,i], type = "p", pch = 19, col = "blue")}
    if (all(i>1, ShowPrevious)){
      if (ShowInputs){graphics::points(xaxis, Inputs[,i-1], type = "p", col = "red")}
      if (ShowOutputs){graphics::points(xaxis, Outputs[,i-1], type = "p", col = "blue")}
    }

    if (secondhold > -0.5 & secondhold <= 1e-10){
      readline(prompt=paste0("Iterate:", i, ".     Press [enter] to see next frame"))
    }
    if (secondhold > 1e-10){Sys.sleep(secondhold)}
  }
end


#' A function for plotting the convergence change over a series of iterates.
#'
#' On the x axis is the index of each iterate. On the y axis are different convergence metric values.
#' @export
#' @param Inputs The Inputs matrix returned by FixedPoint.
#' @param Outputs The Outputs matrix returned by FixedPoint.
#' @param Input_Convergence Convergence vector to be used in the figure. Often makes sense for it to be the same as that used in
#' the Fixedpoint function.
#' @param LinesToDraw A vector describing which convergence metrics to draw. This vector can contain anything in c("Sup_Norm", "Euclidean_Norm",
#' "Sum_Of_Absolute_Value_Of_Residuals", "Smallest_Residual", "Input_Convergence"). Here "Input_Convergence" draws the input convergence vector.
#' "Sup_Norm" gives the greatest residual (in absolute value), "Euclidean_Norm" is the eucleadean norm of the residuals, "Sum_Of_Absolute_Value_Of_Residuals"
#' sums up the absolute value of the residuals and "Smallest_Residual" is the smallest residual by absolute value. Clearly if the function takes
#' and returns a scalar (rather than a vector) then all of these are identical.
#' @param LogScale This is whether or not to use a log scale. If so instead of convergence, log(1+convergence) is shown on the y axis. By default this is true.
#' @param FromIterate This describes the first iterate on the figure.
#' @param ToIterate This describes the last iterate on the figure.
#'
#' @return This function returns nothing. It just shows a plot in the console.
#' @examples
#' Inputs = seq(1,10)
#' Function = function(x){ cos(x) }
#' A = FixedPoint(Function, Inputs, Algorithm = "Anderson")
#' ConvergenceFig(A$Inputs, A$Outputs)
#' # And now to only show a median norm:
#' Differences = A$Inputs - A$Outputs
#' Convergence = apply(Differences,2,function(x){median(abs(x))})
#' ConvergenceFig(A$Inputs, A$Outputs, Convergence, LinesToDraw = "Input_Convergence")
"""
A function for plotting the convergence change over a series of iterates. On the x axis is the index of each iterate. On the y axis are different convergence metric values.
### Takes
 * Inputs - The Inputs matrix returned by FixedPoint.
 * Outputs - The Outputs matrix returned by FixedPoint.
 * Input_Convergence - Convergence vector to be used in the figure. Often makes sense for it to be the same as that used in the Fixedpoint function.
 * LinesToDraw - A vector describing which convergence metrics to draw. This vector can contain anything in c("Sup_Norm", "Euclidean_Norm", "Sum_Of_Absolute_Value_Of_Residuals", "Smallest_Residual", "Input_Convergence"). Here "Input_Convergence" draws the input convergence vector. "Sup_Norm" gives the greatest residual (in absolute value), "Euclidean_Norm" is the eucleadean norm of the residuals, "Sum_Of_Absolute_Value_Of_Residuals" sums up the absolute value of the residuals and "Smallest_Residual" is the smallest residual by absolute value. Clearly if the function takes and returns a scalar (rather than a vector) then all of these are identical.
 * LogScale - This is whether or not to use a log scale. If so instead of convergence, log(1+convergence) is shown on the y axis. By default this is true.
 * FromIterate - This describes the first iterate on the figure.
 * ToIterate - This describes the last iterate on the figure.
### Returns
 * This function returns nothing. It just shows a plot in the console.
### Examples
Inputs = seq(1,10)
Function = function(x){ cos(x) }
A = FixedPoint(Function, Inputs, Algorithm = "Anderson")
ConvergenceFig(A$Inputs, A$Outputs)
 # And now to only show a median norm:
Differences = A$Inputs - A$Outputs
Convergence = apply(Differences,2,function(x){median(abs(x))})
ConvergenceFig(A$Inputs, A$Outputs, Convergence, LinesToDraw = "Input_Convergence")
"""
function ConvergenceFig(Inputs, Outputs, Input_Convergence = c(), LinesToDraw = c("Sup_Norm", "Euclidean_Norm","Sum_Of_Absolute_Value_Of_Residuals", "Smallest_Residual"), LogScale = TRUE, FromIterate = 1, ToIterate = c()){
  if (!is.null(Input_Convergence)){
    if (sum(size(Inputs) != size(Outputs)) > 0.5 | size(Inputs)[2] != length(Input_Convergence)) {
      stop("The matrix of outputs and the matrix of inputs do not have the same shape. They must also have the same width as the length of the ConvergenceVector. As there are some differences in this case nothing can be drawn.")}
  } else {
    if (sum(size(Inputs) != size(Outputs)) > 0.5) {
      stop("The matrix of outputs and the matrix of inputs do not have the same shape. They must also have the same width as the length of the ConvergenceVector. As there are some differences in this case nothing can be drawn.")}
  }

  Lines = list()
  if (is.null(ToIterate)){ToIterate = size(Inputs)[2]}
  xaxis = FromIterate:ToIterate
  Inputs  = matrix(Inputs[, xaxis], ncol = length(xaxis))
  Outputs = matrix(Outputs[, xaxis], ncol = length(xaxis))
  if (!is.null(Input_Convergence)){
    LinesToDraw = unique(c(LinesToDraw, "Input_Convergence"))
    Input_Convergence = Input_Convergence[xaxis]
  }

  if (sum(LinesToDraw != "ConvergenceVector")>0.5){
    ResidualVector = Outputs - Inputs
    Lines$Sup_Norm = apply(ResidualVector , 2, function(x) {max(abs(x))}   )
    Lines$Euclidean_Norm = apply(ResidualVector , 2, function(x) {sqrt(sum((x)^2))}   )
    Lines$Sum_Of_Absolute_Value_Of_Residuals = apply(ResidualVector , 2, function(x) {sum(abs(x))}   )
    Lines$Smallest_Residual = apply(ResidualVector , 2, function(x) {min(abs(x))}   )
    Lines$Input_Convergence = Input_Convergence
    Lines = Lines[LinesToDraw]
  }
  Unlisting = unlist(Lines)
  ytop = max(Unlisting)
  if (LogScale){ybot = 1e-16} else {ybot = 0}


  colours = c("black", "blue", "green", "red", "brown")
  Subtitle = paste0(LinesToDraw[1]," in ", colours[1])
  if ( length(LinesToDraw) > 1.5){
    Subtitle = paste0(Subtitle, ", ")
    for (i in 2:length(LinesToDraw)){
      Subtitle = paste0(Subtitle, LinesToDraw[i]," in ", colours[i])
      if (i %% 2 == 1 & i < length(LinesToDraw)){Subtitle = paste0(Subtitle,", ")}
      if (i %% 2 == 0 & i < length(LinesToDraw)){Subtitle = paste0(Subtitle, "," , "\n")}
    }
  }

  if (LogScale){
    graphics::plot(c(min(xaxis), max(xaxis)), c(ytop, ybot), type = "p", col = 0, xlab = "", ylab = "log(1+Convergence)", main = "Convergence per Iteration", sub = Subtitle, log = "y")
  } else {
    graphics::plot(c(min(xaxis), max(xaxis)), c(ytop, ybot), type = "p", col = 0, xlab = "", ylab = "Convergence", main = "Convergence per Iteration", sub = Subtitle)
  }
  for (line in 1:length(LinesToDraw)){
    if (LogScale){
      graphics::lines(xaxis, log(1+Lines[LinesToDraw[line]][[1]]), type = "l", col = colours[line] )
      graphics::points(xaxis, log(1+Lines[LinesToDraw[line]][[1]]), type = "p", col = colours[line] )
    } else {
      graphics::lines(xaxis, Lines[LinesToDraw[line]][[1]], type = "l", col = colours[line] )
      graphics::points(xaxis, Lines[LinesToDraw[line]][[1]], type = "p", col = colours[line] )
    }
  }
end

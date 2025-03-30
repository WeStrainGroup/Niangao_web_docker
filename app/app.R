library(shiny) 
library(shinyjs) # modify the interface
library(DT) # build datatables
library(zip) # compress files into zip
library(sangerseqR) # ab1 files processing
library(RcppRoll) # sliding windows
library(seqinr) # sequence processing
library(sys) # for executing system commands

# UI
ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    htmlTemplate("www/index.html"),
  ),
  
  tags$style("
      .progress-bar {
        height: 10px;
        font-size: 10px;
        background-color: #8c8e93;
        line-height: 1;
      }
      .shiny-file-input-progress {
        height: 10px;
        border-radius: 0 0 5px 5px;
      }
      .form-control {
        background-color: white !important;
      }
  ")
)

# Server
server <- function(input, output, session) {
  # Generate unique ID and work dir for each session
  session_id <- session$token
  session_dir <- file.path(tempdir(), session_id)
  dir.create(session_dir)

  # Delete the whole dir when session ends
  session$onSessionEnded(function() {
    unlink(session_dir, recursive = TRUE)
  })

  # insert elements into interface
  observe({
    insertUI(
      selector = "#file-input-container",
      ui = fileInput("file", 
                     label=NULL,
                     multiple = TRUE,
                     accept = (".ab1"),
                     width = 1100,
      )
    )
    
    insertUI(
      selector = "#download-button",
      ui = downloadButton("download_result", NULL, style = "border:none; background-color:transparent")
    )
    
    insertUI(
      selector = "#BLAST-or-not",
      ui = selectInput(
        inputId = paste0("BLAST-or-not", input$add_select),
        label = "Do BLAST?",
        choices = c("YES", "NO"),
        width = 300,
      )
    )
    
    insertUI(
      selector = "#database-choosing",
      ui = selectInput(
        inputId = paste0("database-choosing", input$add_select),
        label = "Which database to align?",
        choices = c("16S_ribosomal_RNA", "ITS_RefSeq_Fungi", "ITS_eukaryote_sequences", "28S_fungal_sequences", "18S_fungal_sequences"),
        width = 400,
      )
    )
    
    insertUI(
      selector = "#hits-amount",
      ui = selectInput(
        inputId = paste0("hits-amount", input$add_select),
        label = "How many hits to present?",
        choices = c("1", "5", "10", "20"),
        width = 400,
      )
    )
  })
  
  # define some variables and containers
  ab1_count <- 0 # counting ab1 files
  fasta_all <- character() # prepare to present (fasta sequences integrated)
  processing_notes <- character() # prepare to present (processing notes)
  result_clean <- data.frame(file_name = character(), 
                             length_raw = numeric(), 
                             length_clean = numeric(), 
                             average_quality_raw = numeric(),
                             average_quality_clean = numeric(),
                             deleted_count = numeric(),
                             deleted_ratio = numeric(),
                             deleted_count1 = numeric(),
                             deleted_count2 = numeric(),
                             deleted_count3 = numeric())
  
  # get choices from the interface
  values <- reactiveValues(
    blast_or_not = "YES",
    database_choosing = "16S_ribosomal_RNA",
    hits_amount = "1"
  )

  observe({
    values$blast_or_not <- input[[paste0("BLAST-or-not", input$add_select)]]
    values$database_choosing <- input[[paste0("database-choosing", input$add_select)]]
    values$hits_amount <- input[[paste0("hits-amount", input$add_select)]]
  })

  # set directories of databases
  ITS_RefSeq_Fungi_dir <- "/opt/conda/envs/myenv/db/ITS_RefSeq_Fungi"
  ITS_eukaryote_sequences_dir <- "/opt/conda/envs/myenv/db/ITS_eukaryote_sequences"
  I6S_ribosomal_RNA_dir <- "/opt/conda/envs/myenv/db/16S_ribosomal_RNA"
  Z8S_fungal_sequences_dir <- "/opt/conda/envs/myenv/db/28S_fungal_sequences"
  I8S_fungal_sequences_dir <- "/opt/conda/envs/myenv/db/18S_fungal_sequences"

  # when input is ready, starts processing
  observeEvent(input$file, {
    req(input$file) # ensure the input is not empty

    # function: reset file folders
    reset_dir <- function(dir_path) {
      if (dir.exists(dir_path)) {
        file.remove(list.files(dir_path, full.names = TRUE))
      } else {
        dir.create(dir_path, recursive = TRUE)
      }
    }
    
    # create and reset file folders
    seq_fwd_single <- file.path(getwd(), "default_fwd_single")
    reset_dir(seq_fwd_single)
    seq_rev_single <- file.path(getwd(), "default_rev_single")
    reset_dir(seq_rev_single)
    seq_nonsuffix_single <- file.path(getwd(), "default_nonsuffix_single")
    reset_dir(seq_nonsuffix_single)
    
    result_niangao <- file.path(session_dir, "result_niangao") # folder that contains all the results
    reset_dir(result_niangao)
    seq_clean_single <- file.path(session_dir, "seq_clean_single") # folder that contain all cleaned fasta files
    reset_dir(seq_clean_single)
    seq_assemble_single <- file.path(session_dir, "seq_assemble_single") # folder that contain all assembled cleaned fasta files
    reset_dir(seq_assemble_single)
    
    # prepare for the failed assembly
    failed_assemblies_name <- character()
    failed_assemblies_fasta <- character()
    
    # process with progress bar
    withProgress(message = "Processing...", max = length(input$file$datapath), {
      
      for (i in 1:length(input$file$datapath)) {
        ab1_files <- input$file$datapath[i]
        ab1_name <- substr(basename(input$file$name[i]), 1, nchar(basename(input$file$name[i])) - 4)
        
        Sys.sleep(0.1)  # add delay to observe the update of progress bar
        incProgress(1, detail = basename(input$file$name[i]))  # update progress bar
        
        # select using suffix and count ab1 files
        ext <- tools::file_ext(input$file$name[i])
        if (ext == "ab1") {
          ab1_count <- ab1_count + 1
          
          # extract infomation from ab1s
          abif_data <- sangerseqR::read.abif(ab1_files)
          quality_values <- abif_data@data$PCON.1
          base_calls <- abif_data@data$PBAS.1
          base_calls_split <- strsplit(base_calls, "")[[1]]
          if (length(quality_values) == length(base_calls_split)) {
            paired_data <- data.frame(base_call = base_calls_split, quality_value = quality_values) 
          } else {
            quality_values <- abif_data@data$PCON.2
            base_calls <- abif_data@data$PBAS.2
            base_calls_split <- strsplit(base_calls, "")[[1]]
            paired_data <- data.frame(base_call = base_calls_split, quality_value = quality_values) 
          }
          
          # start clean process
          # clean 1: use sliding windows to filter average quality values
          window_size <- 20 # an appropriate parameter for window size
          quality_averages <- roll_mean(quality_values, n = window_size, align = "center", fill = 0)
          
          filter_size <- 45 # an appropriate parameter for quality value threshold
          windows <- data.frame(num = seq_along(quality_averages), quality_average = quality_averages - filter_size) 
          windows
          
          # add small gaps back to the qualified sequence (gap length <= 10 bases)
          transition_indices <- which(diff(sign(windows$quality_average)) != 0) + 1 # find the positions of sign change after the average quality value reduced by 45
          transition_indices
          
          if (length(transition_indices) > 2) { 
            for (i in seq(2, length(transition_indices) - 1, by = 2)) { # find the positions ofs sign change from positive to negative
              positive_to_negative_position <- windows$num[transition_indices[i]]
              negative_to_positive_position <- windows$num[transition_indices[i+1]]
              position_diff <- negative_to_positive_position - positive_to_negative_position # calculate the length of every gap
              if (position_diff <= 10) { # if the length <= 10, give it 0 value to pass the windows_filter
                windows$quality_average[positive_to_negative_position:negative_to_positive_position] <- 0
              } 
            }
          } 
          
          windows_filter <- subset(windows, quality_average >= 0)
          windows_filter
          
          # if all bases are filtered out, skip other clean steps
          if (nrow(windows_filter) == 0) {
            deleted_count1 = nrow(windows)
            deleted_count2 = 0
            deleted_count3 = 0
            paired_data_clean3 = data.frame(base_calls = character(), quality_values = numeric())
          } else {
            deleted_count1 <- nrow(windows) - (nrow(windows_filter) + window_size - 1) # calculate the number of bases deleted by clean 1
            
            # clean 2 find the longest continuous partition of entire sequence
            x <- windows_filter$num
            current_num <- vector("integer")
            continuous_num <- vector("integer")
            
            for (i in 1:length(x)) {
              if (i > 1 && x[i] == x[i - 1] + 1) {
                current_num <- c(current_num, x[i])
              } else {
                if (length(current_num) > length(continuous_num)) {
                  continuous_num <- current_num
                }
                current_num <- c(x[i]) # move to the next continuous partition
              }
            }
            
            if (length(current_num) > length(continuous_num)) { # in case the last partition is the longest
              continuous_num <- current_num
            }
            
            paired_data_clean2 <- paired_data[(min(continuous_num) - (window_size/2)):(max(continuous_num) + (window_size/2)), ]
            deleted_count2 <- length(windows_filter$num) - length(continuous_num) # calculate the number of bases deleted by clean 2
            
            # clean 3 filter out bases with quality values under 10
            threshold <- 10
            paired_data_clean3 <- paired_data_clean2[paired_data_clean2$quality_value >= threshold, ]
            deleted_count3 <- nrow(paired_data_clean2) - nrow(paired_data_clean3) # calculate the number of bases deleted by clean 3
            
            # generate process report in table form
            result_row <- data.frame(file_name = ab1_name, 
                                     len._raw = nrow(paired_data),
                                     len._trimmed = nrow(paired_data_clean3),
                                     ave._quality_raw = round(mean(paired_data$quality_value), 1), 
                                     ave._quality_trimmed = round(mean(paired_data_clean3$quality_value), 1),
                                     len._deleted = nrow(paired_data) - nrow(paired_data_clean3),
                                     ratio_deleted = round((nrow(paired_data) - nrow(paired_data_clean3))/nrow(paired_data) * 100, 1),
                                     deleted_count1 = deleted_count1,
                                     deleted_count2 = deleted_count2,
                                     deleted_count3 = deleted_count3)
            result_clean <- rbind(result_clean, result_row)
            
            # generate fasta
            fasta_single <- c(">", ab1_name, "\n", paste0(paired_data_clean3$base_call, collapse = ""))
            fasta_all <- c(fasta_all, fasta_single, "\n")
            
            write.fasta(paired_data_clean3$base_call, ab1_name, file.path(seq_clean_single, paste0(ab1_name, ".fasta")))
          }
        }
      }
      
      # prepare for CAP3 assembly
      consensusSeqs <- reactiveVal(list()) # consensus sequences from CAP3 assembly
      assemble_all <- character() # successfully assembled fasta
      ace_content_all <- character() # ace file from CAP3 assembly
      fwd_only_all <- character() # only fwd file is input in a pair
      rev_only_all <- character() # only rev file is input in a pair
      nonsuffix_all <- character() # nonsuffix file is input while paired fasta exists
      
      # separate paired fasta and nonsuffix fasta
      fasta_index <- grep("^>", fasta_all)
      fasta_headers <- fasta_all[sort(fasta_index + 1)]
      
      fwd_headers <- grep("_Fwd", fasta_headers, value = TRUE)
      rev_headers <- grep("_Rev", fasta_headers, value = TRUE)
      nonsuffix_headers <- setdiff(fasta_headers, c(fwd_headers, rev_headers))
      
      fwd_samples <- gsub("^>|_Fwd.*", "", fwd_headers)
      rev_samples <- gsub("^>|_Rev.*", "", rev_headers)
      
      paired_samples <- union(fwd_samples, rev_samples)
      
      # if no paired fasta, skip assembly process
      if (length(paired_samples) == 0) {
        # set a parameter (used in the output session)
        single_or_paired <- "single"
      } else {
        # set a parameter (used in the output session)
        single_or_paired <- "paired"
        
        # save nonsuffix samples into folder first
        for (header in nonsuffix_headers) {
          header_index <- which(fasta_all == header)
          if (length(header_index) > 0) {
            seq <- fasta_all[header_index + 2]
            sample_name_nonsuffix <- gsub("^>", "", header)
            nonsuffix_fasta_single <- c(">", sample_name_nonsuffix, "\n", paste0(seq), "\n")
            nonsuffix_all <- c(nonsuffix_all, nonsuffix_fasta_single)
            
            # write nonsuffix fasta
            write.fasta(strsplit(seq, "")[[1]], sample_name_nonsuffix, file.path(seq_nonsuffix_single, paste0(sample_name_nonsuffix, ".fasta")))
          }
        }
        
        # function: get header and seq for a sample
        fasta_pairs <- lapply(paired_samples, function(sample) {
          
          fwd_header <- ""
          rev_header <- ""
          fwd_seq <- ""
          rev_seq <- ""
          
          fwd_header <- fwd_headers[grep(sample, fwd_headers)]
          rev_header <- rev_headers[grep(sample, rev_headers)]
          
          fwd_seq_index <- which(fasta_all == fwd_header)
          rev_seq_index <- which(fasta_all == rev_header)
          
          if (length(fwd_header) != 0){
            fwd_seq <- fasta_all[fwd_seq_index + 2]
          } 
          if (length(rev_header) != 0) {
            rev_seq <- fasta_all[rev_seq_index + 2]
          }
          
          return(list(sample = sample, fwd = fwd_seq, rev = rev_seq, fwd_header = fwd_header, rev_header = rev_header))
        })
        
        # process paired samples
        sample_names <- character()
        
        for (pair in fasta_pairs) {
          fwd_seq <- pair$fwd
          rev_seq <- pair$rev
          fwd_header <- pair$fwd_header
          rev_header <- pair$rev_header
          sample_name <- pair$sample
          sample_names <- c(sample_names, sample_name)
          mergedFile <- tempfile(fileext = ".fasta") # merge paired fasta to use CAP3
          outFile <- tempfile(fileext = ".ace") # temporarily save output ace of CAP3
          
          if (length(fwd_header) == 0) {
            rev_only_single <- c(">", sample_name,"_Rev\n", paste0(rev_seq), "\n", collapse = "")
            rev_only_all <- c(rev_only_all, rev_only_single)
            
            # write rev_single_fasta
            write.fasta(sequences = strsplit(rev_seq, "")[[1]], 
                        names = paste0(sample_name, "_Rev"), 
                        file.out = file.path(seq_rev_single, paste0(sample_name, "_Rev", ".fasta")))
            
          } else if (length(rev_header) == 0) {
            fwd_only_single <- c(">", sample_name,"_Fwd\n", paste0(fwd_seq), "\n", collapse = "")
            fwd_only_all <- c(fwd_only_all, fwd_only_single)
            
            # write fwd_single_fasta
            write.fasta(sequences = strsplit(fwd_seq, "")[[1]],
                        names = paste0(sample_name, "_Fwd"),
                        file.out = file.path(seq_fwd_single, paste0(sample_name, "_Fwd", ".fasta")))
            
          } else {
            # find CAP3.exe
            cap3_path <- "./cap3"
            
            # merge paired fasta to input for CAP3
            merge_seq <- paste(">", sample_name, "_Fwd\n", fwd_seq, "\n",
                               ">", sample_name, "_Rev\n", rev_seq, sep="")
            writeLines(merge_seq, mergedFile)
            
            # run CAP3
            system2(cap3_path, 
                    args = mergedFile, 
                    stdout = outFile, 
                    stderr = "output_log.txt")
            
            # function: extract consensus sequence
            extractConsensus <- function(ace_content) {
              consensus <- list()
              current_contig <- NULL
              seq_lines <- c()
              
              for (line in ace_content) {
                if (grepl("Contig", line)) {
                  if (!is.null(current_contig) && length(seq_lines) > 0) {
                    consensus[[current_contig]] <- paste(seq_lines, collapse = "")
                  }
                  current_contig <- sub(".*Contig (\\d+).*", "\\1", line)
                  seq_lines <- c()
                } else if (grepl("consensus", line)) {
                  seq_lines <- c(seq_lines, gsub("consensus\\s+", "", line))
                }
              }
              
              if (!is.null(current_contig) && length(seq_lines) > 0) {
                consensus[[current_contig]] <- paste(seq_lines, collapse = "")
              } 
              
              return(list(consensus = consensus, ace_content = ace_content))
            }
            
            # extract ace file content in CAP3 results
            ace_content <- readLines(outFile)
            cap3_results <- extractConsensus(ace_content)
            
            # check whether succeed in assembly
            if (length(cap3_results$consensus) == 0) {
              failed_assemblies_name <- c(failed_assemblies_name, sample_name) # Record failed sample
              
              # select the longer sequence in a pair to output
              if (nchar(fwd_seq) >= nchar(rev_seq)) {
                longer_seq <- fwd_seq
                seq_type <- "Fwd" 
              } else {
                longer_seq <- rev_seq
                seq_type <- "Rev"
              } 
              
              # write longer seq into fasta
              write.fasta(longer_seq, paste0(sample_name, "_", seq_type), file.path(seq_assemble_single, paste0(sample_name, "_", seq_type, ".fasta")))
              failed_assemblies_fasta <- c(failed_assemblies_fasta,
                                           paste0(">", sample_name, "_", seq_type, "\n", longer_seq, "\n"))
            } else {
              # update current_consensus 
              current_consensus <- consensusSeqs()
              current_consensus[[sample_name]] <- cap3_results$consensus
              
              # save all ace file content
              ace_content_all <- c(ace_content_all, cap3_results$ace_content)
              
              # generate assembled fasta
              assemble_single <- c(">", sample_name, "_assemble", "\n", paste0(unlist(cap3_results$consensus), "\n", collapse = ""))
              assemble_all <- c(assemble_all, assemble_single)
              
              # write assembled fasta
              write.fasta(current_consensus[[sample_name]], sample_name, file.path(seq_assemble_single, paste0(sample_name, "_assemble", ".fasta")))
              
              # update consensus seq to interface
              consensusSeqs(current_consensus)
            }
          }
        }
      }

      # output to interface
      # processing_notes
      if (ab1_count == length(input$file$datapath)) {
        processing_notes <- "All files were successfully processed!<br> You can download the results using the given button."
      } else {
        processing_notes <- paste0("The uploaded files contain invalid file types, and only *.ab1 files were processed.")
      }
      
      # Process failed assemblies notification
      if (length(failed_assemblies_name) > 0) {
        processing_notes <- paste0(processing_notes, 
                                   "<br><strong>Assembly failed for the following samples:</strong><br>",
                                   paste(failed_assemblies_name, collapse = ", "), "<br>")
      }
      
      # display in notebox
      shinyjs::html("notebox", processing_notes)
      shinyjs::show("notebox")
      
      # show blast is running
      if (values$blast_or_not == "YES") {
        blast_notes <- "<br><strong>BLAST is running...</strong>"
        shinyjs::html("BLAST", blast_notes)
      }

      # trimming report
      observe({
        shinyjs::html("report", "")
        insertUI(
          selector = "#report",
          ui = dataTableOutput("result_clean")
        )
      })
      output$result_clean <- DT::renderDataTable({
        DT::datatable(result_clean[, 1:7],
                      colnames = c("File Name", "Len. Raw", "Len. Trimmed", "Ave. Quality Raw", 
                                 "Ave. Quality Trimmed", "Deleted Bases", "Deleted Ratio (%)")
        )
      })
      
      # single or paired?
      if (single_or_paired == "single") {
        # show cleaned fasta on the interface
        observe({
          shinyjs::html("fasta-clean", paste0(fasta_all, collapse=""))
        })
      } else {
        # show cleaned fasta on the interface
        observe({
          consensus <- consensusSeqs()
          # get every consensus sequence
          output_text <- lapply(names(consensus), function(sample_name) {
            paste0(">", sample_name, "_assemble", "<br/>", paste(unlist(consensus[[sample_name]]), collapse = ""), "<br/>")
          })
          
          # add rev only, fwd only, failed assembiles and nonsuffixes
          if (length(rev_only_all) > 0) {
            output_text <- c(output_text, rev_only_all)
          }
          if (length(fwd_only_all) > 0) {
            output_text <- c(output_text, fwd_only_all)
          }
          if (length(failed_assemblies_fasta) > 0) {
            output_text <- c(output_text, failed_assemblies_fasta)
          }
          if (length(nonsuffix_all) > 0) {
            output_text <- c(output_text, nonsuffix_all)
          }
          # combine all into onevariable
          output_text_combined <- paste(unlist(output_text), collapse = "")
          
          # show on the interface
          shinyjs::html("fasta-clean", output_text_combined)
        })
        
        # additional downloadable output
        assemble_all <- c(assemble_all, rev_only_all, fwd_only_all, failed_assemblies_fasta, nonsuffix_all)
        writeLines(paste0(assemble_all, collapse = ""), file.path(result_niangao, "seq_trimed_assembled_all.fasta"))
        writeLines(paste0(ace_content_all, collapse = "\n"), file.path(result_niangao, "cap3_aseembly_details.ace"))
        zip(zipfile = file.path(result_niangao, "seq_trimed_assembled_single.zip"),
            files = c(
              list.files(seq_assemble_single, full.names = TRUE, pattern = "*.fasta"),
              list.files(seq_rev_single, full.names = TRUE, pattern = "*.fasta"),
              list.files(seq_fwd_single, full.names = TRUE, pattern = "*.fasta"),
              list.files(seq_nonsuffix_single, full.names = TRUE, pattern = "*.fasta")),
            mode = "cherry-pick")
      }

      # shared downloadable output
      writeLines(paste0(fasta_all, collapse = ""), file.path(result_niangao, "seq_trimed_all.fasta"))
      write.csv(result_clean, file.path(result_niangao, "trimming_report.csv"))
      zip(zipfile = file.path(result_niangao, "seq_trimed_single.zip"), 
          files = list.files(seq_clean_single, full.names = TRUE, pattern = "*.fasta"),
          mode = "cherry-pick")
      
      if (values$blast_or_not == "YES") {
        database_path <- switch(values$database_choosing,
                                    "ITS_RefSeq_Fungi" = ITS_RefSeq_Fungi_dir,
                                    "ITS_eukaryote_sequences" = ITS_eukaryote_sequences_dir,
                                    "16S_ribosomal_RNA" = I6S_ribosomal_RNA_dir,
                                    "28S_fungal_sequences.tar.gz" = Z8S_fungal_sequences_dir,
                                    "18S_fungal_sequences.tar.gz" = I8S_fungal_sequences_dir)  

        # Construct full paths
        query_file <- file.path(result_niangao, if (single_or_paired == "single") "seq_trimed_all.fasta" else "seq_trimed_assembled_all.fasta")
        output_file <- file.path(session_dir, "blast_result.tsv")
        blast_command <- "/opt/conda/envs/myenv/bin/blastn"

        # Get the number of available cores
        num_threads <- parallel::detectCores(logical = TRUE)

        # Run blast
        blast_args <- c("-query", query_file,
                "-db", database_path,
                "-out", output_file,
                sprintf('-outfmt "6 qseqid sseqid staxids pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore"'),
                "-max_target_seqs", values$hits_amount,
                "-num_threads", as.character(num_threads),
                "-evalue", "1e-2")

        blast_result <- system2(blast_command, args = blast_args, stdout = TRUE, stderr = TRUE)

        # Check for errors
        if (!is.null(attr(blast_result, "status")) && attr(blast_result, "status") != 0) {
          error_message <- paste("BLAST failed with error:", paste(blast_result, collapse = "\n"))
          shinyjs::html("BLAST", error_message)  # Display error in the UI
          showModal(modalDialog(
              title = "BLAST Error",
              error_message,
              easyClose = TRUE,
              footer = NULL
          ))
        } else {
          # BLAST result annotation
          lineage_file <- file.path(session_dir, "lineage.txt")
          final_output_file <- file.path(result_niangao, "blast_result_final.tsv")

          # Extract lineage
          taxonkit_cmd <- paste("cut -f3", output_file, "|",
                                "sed 's/;.*//' |",
                                "taxonkit reformat -I 1 -F -P -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' >", lineage_file)
          system2("bash", args = c("-c", shQuote(taxonkit_cmd)))

          # Combine results and add header
          paste_cmd <- paste("paste", output_file, lineage_file, "|",
                            "csvtk add-header -t -n 'qseqid,sseqid,staxids,pident,length,qcovs,mismatch,gapopen,qstart,qend,sstart,send,sstrand,evalue,bitscore,taxids,kingdom,phylum,class,order,family,genus,species' >", final_output_file)
          system2("bash", args = c("-c", shQuote(paste_cmd)))

          # present blast results to the interface
          observe({
            shinyjs::html("BLAST", "")

            insertUI(
              selector = "#BLAST",
              ui = dataTableOutput("BLAST")
            )

            output$BLAST <- DT::renderDataTable({
              blast_result_table <- tryCatch({
                pre_blast_result_table <- read.table(final_output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

                # Present species names without the 's__' prefix
                pre_blast_result_table$species <- sub("^s__", "", pre_blast_result_table$species)

                # Change the column names
                blast_result_table <- pre_blast_result_table[, c("qseqid","species", "pident", "evalue", "bitscore", "length", "qcovs")]
                colnames(blast_result_table) <- c("Input", "Scientific Name", "Per. Ident (%)", "E Value", "Bit Score", "Aligned len.", "Query Cover (%)")
                
                # Present the table
                DT::datatable(blast_result_table)
              }, error = function(e) {
                showModal(modalDialog(
                    title = "File Read Error",
                    paste("Could not read BLAST output file:", e$message),
                    easyClose = TRUE,
                    footer = NULL))
                return(NULL)
              })

              if (!is.null(blast_result_table)) {
                return(blast_result_table)
              }
            })  
          })
        }
      }
      
      # pack to download button
      output$download_result <- downloadHandler(
        filename = function() {
          paste("result_niangao_", format(Sys.Date(), "%Y%m%d"), ".zip", sep = "")
        },
        content = function(file) {
          # Write the README.md file to the result directory
          readme_path <- file.path("www", "README.md")
          readme_content <- readLines(readme_path, warn = FALSE)
          writeLines(readme_content, file.path(result_niangao, "README.md"))

          file.copy(zip(zipfile = file.path(session_dir,"output.zip"), 
                        files = list.files(result_niangao, full.names = TRUE),
                        mode = "cherry-pick"),
                    file)
        }
      )
    })
  })
}

# run the app
options(shiny.host = '0.0.0.0')
options(shiny.port = 3838)

runApp(list(ui = ui, server = server))
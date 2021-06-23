from treasureisland.gi_driver import gi_driver 


driver = gi_driver("C:/Users/USER/GenomicIslandPrediction/treasureisland/genome/test.fasta")

pred = driver.get_predictions()

driver.predictions_to_csv(pred)

driver.predictions_to_excel(pred)

driver.predictions_to_text(pred)


set :application, "ngs-aggregate_results"
set :repo_url, "git@github.com:nebiolabs/ngs-aggregate_results.git"
set :branch, "master"
set :deploy_to, "#{ENV['DEPLOY_BASE']}/ngs-aggregate_results"

role :app, "[#{ENV['DEPLOY_USER']}@#{ENV['DEPLOY_HOST']}]"
append :linked_files, "config/database.yml", "config/bc_purity_database.yml", "config/config.yml", "config/.env", "config/master.key"
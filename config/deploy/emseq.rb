set :application, "EM-seq"
set :repo_url, "git@github.com:nebiolabs/EM-seq.git"
set :branch, ENV["DEPLOY_TAG"]
set :deploy_to, ENV['DEPLOY_BASE']

role :app, ["#{ENV['DEPLOY_USER']}@#{ENV['DEPLOY_HOST']}"]
append :linked_files, "nextflow.config"

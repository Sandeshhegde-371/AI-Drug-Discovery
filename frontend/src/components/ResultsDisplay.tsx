import { Download, CheckCircle, TrendingUp, Award, BarChart3 } from 'lucide-react';
import { Card } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from './ui/table';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ScatterChart, Scatter } from 'recharts';
import type { Pipeline } from './PipelineSelector';

interface ResultsDisplayProps {
  pipeline: Pipeline;
  results: any[];
  fileName: string;
}

export function ResultsDisplay({ pipeline, results, fileName }: ResultsDisplayProps) {
  const downloadResults = () => {
    if (!results.length) return;

    const headers = Object.keys(results[0]);
    const csv = [
      headers.join(','),
      ...results.map(row => headers.map(header => row[header]).join(',')),
    ].join('\n');

    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `results_${fileName}`;
    a.click();
    window.URL.revokeObjectURL(url);
  };

  const getPipelineTitle = () => {
    switch (pipeline) {
      case 'generative':
        return 'Generated Molecules';
      case 'virtual':
        return 'Screening Results';
      case 'properties':
        return 'Property Predictions';
      case 'prioritize':
        return 'Prioritized Candidates';
      default:
        return 'Results';
    }
  };

  const getPropertyDistributionData = () => {
    if (!results.length) return [];
    
    const bins: { [key: string]: number } = {};
    results.forEach(r => {
      const score = Math.floor(r.score * 10) / 10;
      bins[score] = (bins[score] || 0) + 1;
    });

    return Object.entries(bins).map(([score, count]) => ({
      score: parseFloat(score),
      count,
    })).sort((a, b) => a.score - b.score);
  };

  const getScatterData = () => {
    // Find first two numeric columns for scatter plot
    if (!results.length) return [];
    
    const keys = Object.keys(results[0]);
    const numericKeys = keys.filter(key => {
      const val = parseFloat(results[0][key]);
      return !isNaN(val) && key !== 'score';
    }).slice(0, 2);

    return results.slice(0, 50).map(r => ({
      x: parseFloat(r[numericKeys[0]] || 0),
      y: parseFloat(r[numericKeys[1]] || 0),
      score: r.score,
    }));
  };

  const getTableColumns = () => {
    if (!results.length) return [];
    const allKeys = Object.keys(results[0]);
    // Show first 5 columns plus score
    return allKeys.filter(k => k !== 'score').slice(0, 4).concat(['score']);
  };

  const highConfidence = results.filter(r => r.score > 0.8).length;
  const avgScore = (results.reduce((sum, r) => sum + r.score, 0) / results.length).toFixed(3);

  return (
    <Card className="p-8 bg-white/70 backdrop-blur-sm border-purple-100 shadow-lg shadow-purple-100/50">
      <div className="mb-8">
        <div className="flex items-start justify-between mb-6">
          <div className="flex items-center gap-3">
            <div className="w-12 h-12 bg-gradient-to-br from-green-600 to-emerald-600 rounded-lg flex items-center justify-center">
              <CheckCircle className="w-6 h-6 text-white" />
            </div>
            <div>
              <h3 className="text-gray-900">Analysis Complete</h3>
              <p className="text-sm text-gray-600">
                {getPipelineTitle()} • {results.length} results generated
              </p>
            </div>
          </div>
          <Button 
            onClick={downloadResults} 
            className="gap-2 bg-gradient-to-r from-purple-600 to-indigo-600 hover:from-purple-700 hover:to-indigo-700 shadow-lg shadow-purple-500/30"
          >
            <Download className="w-4 h-4" />
            Download CSV
          </Button>
        </div>

        {/* Stats Cards */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-8">
          <div className="bg-gradient-to-br from-blue-50 to-indigo-50 border-2 border-blue-200 rounded-xl p-5">
            <div className="flex items-center gap-3 mb-2">
              <div className="w-10 h-10 bg-blue-600 rounded-lg flex items-center justify-center">
                <BarChart3 className="w-5 h-5 text-white" />
              </div>
              <p className="text-blue-900">Total Results</p>
            </div>
            <p className="text-2xl text-blue-900 ml-13">{results.length}</p>
          </div>
          
          <div className="bg-gradient-to-br from-green-50 to-emerald-50 border-2 border-green-200 rounded-xl p-5">
            <div className="flex items-center gap-3 mb-2">
              <div className="w-10 h-10 bg-green-600 rounded-lg flex items-center justify-center">
                <Award className="w-5 h-5 text-white" />
              </div>
              <p className="text-green-900">High Confidence</p>
            </div>
            <p className="text-2xl text-green-900 ml-13">
              {highConfidence} <span className="text-sm">({((highConfidence / results.length) * 100).toFixed(1)}%)</span>
            </p>
          </div>
          
          <div className="bg-gradient-to-br from-purple-50 to-indigo-50 border-2 border-purple-200 rounded-xl p-5">
            <div className="flex items-center gap-3 mb-2">
              <div className="w-10 h-10 bg-purple-600 rounded-lg flex items-center justify-center">
                <TrendingUp className="w-5 h-5 text-white" />
              </div>
              <p className="text-purple-900">Average Score</p>
            </div>
            <p className="text-2xl text-purple-900 ml-13">{avgScore}</p>
          </div>
        </div>
      </div>

      {/* Charts Section */}
      <div className="space-y-8 mb-8">
        <div className="bg-gradient-to-br from-gray-50 to-slate-50 border border-gray-200 rounded-xl p-6">
          <h4 className="text-gray-900 mb-4 flex items-center gap-2">
            <BarChart3 className="w-5 h-5 text-purple-600" />
            Score Distribution
          </h4>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={getPropertyDistributionData()}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
              <XAxis dataKey="score" stroke="#6b7280" />
              <YAxis stroke="#6b7280" />
              <Tooltip 
                contentStyle={{ 
                  backgroundColor: 'white', 
                  border: '1px solid #e5e7eb',
                  borderRadius: '8px',
                  boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)'
                }} 
              />
              <Legend />
              <Bar dataKey="count" fill="url(#colorGradient)" radius={[8, 8, 0, 0]} />
              <defs>
                <linearGradient id="colorGradient" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="0%" stopColor="#9333ea" />
                  <stop offset="100%" stopColor="#6366f1" />
                </linearGradient>
              </defs>
            </BarChart>
          </ResponsiveContainer>
        </div>

        {getScatterData().length > 0 && (
          <div className="bg-gradient-to-br from-gray-50 to-slate-50 border border-gray-200 rounded-xl p-6">
            <h4 className="text-gray-900 mb-4 flex items-center gap-2">
              <TrendingUp className="w-5 h-5 text-purple-600" />
              Data Properties (Top 50)
            </h4>
            <ResponsiveContainer width="100%" height={280}>
              <ScatterChart>
                <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
                <XAxis dataKey="x" name="Property 1" stroke="#6b7280" />
                <YAxis dataKey="y" name="Property 2" stroke="#6b7280" />
                <Tooltip 
                  cursor={{ strokeDasharray: '3 3' }}
                  contentStyle={{ 
                    backgroundColor: 'white', 
                    border: '1px solid #e5e7eb',
                    borderRadius: '8px',
                    boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)'
                  }}
                />
                <Legend />
                <Scatter name="Data Points" data={getScatterData()} fill="#8b5cf6" />
              </ScatterChart>
            </ResponsiveContainer>
          </div>
        )}
      </div>

      {/* Results Table */}
      <div className="bg-gradient-to-br from-gray-50 to-slate-50 border border-gray-200 rounded-xl overflow-hidden">
        <div className="p-6 pb-0">
          <h4 className="text-gray-900 mb-4 flex items-center gap-2">
            <Award className="w-5 h-5 text-purple-600" />
            Top 10 Results
          </h4>
        </div>
        <div className="overflow-x-auto">
          <Table>
            <TableHeader>
              <TableRow className="bg-white/50">
                <TableHead className="font-semibold">Rank</TableHead>
                {getTableColumns().map((col, idx) => (
                  <TableHead key={idx} className="capitalize font-semibold">
                    {col.replace(/_/g, ' ')}
                  </TableHead>
                ))}
                <TableHead className="font-semibold">Status</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {results.slice(0, 10).map((result, index) => (
                <TableRow key={index} className="hover:bg-white/50 transition-colors">
                  <TableCell>
                    <div className="flex items-center gap-2">
                      <div className={`w-8 h-8 rounded-lg flex items-center justify-center text-sm ${
                        index < 3 ? 'bg-gradient-to-br from-purple-600 to-indigo-600 text-white' : 'bg-gray-100 text-gray-600'
                      }`}>
                        #{index + 1}
                      </div>
                    </div>
                  </TableCell>
                  {getTableColumns().map((col, idx) => (
                    <TableCell key={idx} className="max-w-[200px] truncate">
                      {col === 'score' ? (
                        <span className={`font-semibold ${result[col] > 0.8 ? 'text-green-600' : result[col] > 0.6 ? 'text-blue-600' : 'text-gray-600'}`}>
                          {typeof result[col] === 'number' ? result[col].toFixed(3) : result[col]}
                        </span>
                      ) : (
                        <span className="text-gray-700">{typeof result[col] === 'number' ? parseFloat(result[col]).toFixed(2) : result[col]}</span>
                      )}
                    </TableCell>
                  ))}
                  <TableCell>
                    {result.score > 0.8 ? (
                      <Badge className="bg-gradient-to-r from-green-600 to-emerald-600 text-white border-0">High</Badge>
                    ) : result.score > 0.6 ? (
                      <Badge variant="secondary" className="bg-blue-100 text-blue-700 border-blue-200">Medium</Badge>
                    ) : (
                      <Badge variant="outline" className="border-gray-300 text-gray-600">Low</Badge>
                    )}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </div>
      </div>
    </Card>
  );
}
